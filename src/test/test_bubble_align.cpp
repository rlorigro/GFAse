#include "IncrementalIdMap.hpp"
#include "gfa_to_handle.hpp"
#include "graph_utility.hpp"
#include "ContactGraph.hpp"
#include "Bipartition.hpp"
#include "hash_graph.hpp"
#include "Filesystem.hpp"
#include "Sequence.hpp"
#include "Hasher2.hpp"
#include "Timer.hpp"
#include "CLI11.hpp"
#include "Sam.hpp"
#include "Bam.hpp"
#include "minimap.h"

#include "SvgPlot.hpp"

using gfase::for_element_in_sam_file;
using gfase::unpaired_mappings_t;
using gfase::paired_mappings_t;
using gfase::contact_map_t;
using gfase::AlignmentBlock;
using gfase::AlignmentChain;

using gfase::gfa_to_handle_graph;
using gfase::IncrementalIdMap;
using gfase::ContactGraph;
using gfase::Node;
using gfase::Bipartition;
using gfase::SamElement;
using gfase::Sequence;
using gfase::Hasher2;
using gfase::Bubble;
using gfase::Timer;
using gfase::Bam;

using bdsg::HashGraph;
using ghc::filesystem::path;
using CLI::App;

#include <unordered_map>

using std::unordered_map;


void print_minimap_alignment_block(mm_mapopt_t& map_options, mm_idx_t* mi, mm_reg1_t* r2, const string& name, const string& query){
    string type;
    if (r2->id == r2->parent) type = r2->inv? 'I' : 'P';
    else type = r2->inv? 'i' : 'S';

    assert(r2->p); // with MM_F_CIGAR, this should not be NULL
    printf("%s\t%lu\t%d\t%d\t%c\t", name.c_str(), query.size(), r2->qs, r2->qe, "+-"[r2->rev]);
    printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\tTP:%s\tcg:Z:", mi->seq[r2->rid].name, mi->seq[r2->rid].len, r2->rs, r2->re, r2->mlen, r2->blen, r2->mapq, type.c_str());

    for (size_t i = 0; i < r2->p->n_cigar; ++i) // IMPORTANT: this gives the CIGAR in the aligned regions. NO soft/hard clippings!
        printf("%d%c", r2->p->cigar[i] >> 4, MM_CIGAR_STR[r2->p->cigar[i] & 0xf]);
    putchar('\n');
}


AlignmentChain map_sequences(
        const vector<string>& targets,
        const vector<string>& target_names,
        const vector<string>& queries,
        const vector<string>& query_names){

    AlignmentChain result;

    vector<const char*> c_targets;
    vector<const char*> c_names;

    c_targets.reserve(targets.size());
    c_names.reserve(targets.size());

    for(size_t i=0; i<targets.size(); i++) {
        c_targets.push_back(targets[i].c_str());
        c_names.push_back(target_names[i].c_str());
    }

    mm_idxopt_t index_options;
    mm_mapopt_t map_options;

    mm_verbose = 3; // disable message output to stderr
    mm_set_opt(0, &index_options, &map_options);
    mm_set_opt("asm10", &index_options, &map_options);

    index_options.k = 21;
    map_options.flag |= MM_F_CIGAR; // perform alignment
    map_options.flag |= MM_F_EQX;

    mm_idx_t *mi = mm_idx_str(
            index_options.w,
            index_options.k,
            int(0),
            index_options.bucket_bits,
            int(targets.size()),
            c_targets.data(),
            c_names.data()
    );

    mm_tbuf_t *tbuf = mm_tbuf_init(); // thread buffer; for multi-threading, allocate one tbuf for each thread
    for (size_t i=0; i<queries.size(); i++) {
        auto& query = queries[i];
        auto& name = query_names[i];

//        cerr << "ALIGNING" << '\n';
//        cerr << target_names[0] << " " << query_names[i] << ' ' << query.size() << ' ' << targets[0].size() << '\n';

        mm_mapopt_update(&map_options, mi); // this sets the maximum minimizer occurrence; TODO: set a better default in mm_mapopt_init()!

        int n_reg;
        mm_reg1_t *reg;
        reg = mm_map(mi, query.size(), query.c_str(), &n_reg, tbuf, &map_options, name.c_str()); // get all hits for the query

        for (int j = 0; j < n_reg; ++j) { // traverse hits
            mm_reg1_t *r2 = &reg[j];

            assert(r2->p); // with MM_F_CIGAR, this should not be NULL

            if (r2->id == r2->parent){
                AlignmentBlock block(
                        r2->rs,
                        r2->re,
                        r2->qs,
                        r2->qe,
                        0,
                        0,
                        0,
                        0,
                        r2->rev);

                for (uint32_t k = 0; k < r2->p->n_cigar; ++k) { // IMPORTANT: this gives the CIGAR in the aligned regions. NO soft/hard clippings!
                    uint32_t length = r2->p->cigar[k] >> 4;
                    char operation = MM_CIGAR_STR[r2->p->cigar[k] & 0xf];

                    if (operation == '='){
                        block.n_matches += length;
                    }
                    else if (operation == 'X'){
                        block.n_mismatches += length;
                    }
                    else if (operation == 'I'){
                        block.n_inserts += length;
                    }
                    else if (operation == 'D'){
                        block.n_deletes += length;
                    }
                }

                result.chain.emplace_back(block);

                free(r2->p);
            }
        }
        free(reg);
    }
    mm_tbuf_destroy(tbuf);
    mm_idx_destroy(mi);

    result.sort_chains(true);

    return result;
}


void phase_hic(path output_dir, path gfa_path, size_t n_threads){
    Timer t;

    if (exists(output_dir)){
        throw runtime_error("ERROR: output directory exists already");
    }
    else {
        create_directories(output_dir);
    }

//    write_config(output_dir, sam_path, required_prefix, min_mapq, n_threads);

    // Id-to-name bimap for reference contigs
    IncrementalIdMap<string> id_map(false);

    GfaReader reader(gfa_path);

    // TODO: move this into domain of bdsg graph instead of GFA reader
    vector<Sequence> sequences;
    unordered_map<string,size_t> name_to_sequence;
    reader.for_each_sequence([&](string& name, string& sequence){
        id_map.try_insert(name);
        name_to_sequence.emplace(name, sequences.size());
        sequences.emplace_back(name, sequence);
    });

    double sample_rate = 0.04;
    size_t k = 22;
    size_t n_iterations = 6;

    Hasher2 hasher(k, sample_rate, n_iterations, n_threads);

    hasher.hash(sequences);
    hasher.write_results(output_dir);

    size_t max_hits = 5;
    double min_similarity = 0.2;

    unordered_map <pair <string,string>, AlignmentChain> alignments;
    unordered_set <pair <string,string> > visited_pairs;

    ContactGraph alignment_graph;

    hasher.for_each_overlap(max_hits, min_similarity,[&](const string& a, const string& b, int64_t n_hashes, int64_t total_hashes){
        // Skip self hits
        if (a == b){
            return;
        }

        auto& seq_a = sequences[name_to_sequence[a]].sequence;
        auto& seq_b = sequences[name_to_sequence[b]].sequence;

        pair<string,string> ordered_pair;
        if (seq_a.size() > seq_b.size()){
            ordered_pair = {a,b};
        }
        else{
            ordered_pair = {b,a};
        }

        auto result = visited_pairs.find(ordered_pair);

        if (result == visited_pairs.end()) {
            visited_pairs.emplace(ordered_pair);
        }
        else{
            // Don't align pairs twice
            return;
        }

        // Longer length is first
        auto length_a = int32_t(max(seq_a.size(), seq_b.size()));
        auto length_b = int32_t(min(seq_a.size(), seq_b.size()));

        double size_ratio = double(length_b) / double(length_a);

        if (size_ratio < min_similarity){
            // Don't align reads with a size_ratio that would make min_similarity impossible during alignment
            // Occasionally needed where hash similarity fails due to repetitiveness
            return;
        }

        AlignmentChain alignment;

        if (seq_a.size() > seq_b.size()) {
            const vector<string>& targets = {seq_a};
            const vector<string>& target_names = {a};
            const vector<string>& queries = {seq_b};
            const vector<string>& query_names = {b};

            alignment = map_sequences(targets, target_names, queries, query_names);
        }
        else{
            const vector<string>& targets = {seq_b};
            const vector<string>& target_names = {b};
            const vector<string>& queries = {seq_a};
            const vector<string>& query_names = {a};

            alignment = map_sequences(targets, target_names, queries, query_names);
        }

        if (not alignment.empty()) {
            cerr << ordered_pair.first << ' ' << ordered_pair.second << '\n';
            for (auto& item: alignment.chain) {
                cerr << item.get_reversal_char() << '\t' << '(' << item.ref_start << ',' << item.ref_stop << ")\t("
                     << item.query_start << ',' << item.query_stop << ')' << '\t' << item.get_identity() << '\n';
            }
            cerr << '\n';

            // Make sure to retain the ordering by size
            auto id_a = int32_t(id_map.get_id(ordered_pair.first));
            auto id_b = int32_t(id_map.get_id(ordered_pair.second));

            // Clip maximum matches to the length of the longer node
            auto total_matches = min(length_a, int32_t(alignment.get_approximate_non_overlapping_matches()));

            auto alignment_coverage = double(total_matches) / double(length_a);

            if (alignment_coverage < min_similarity){
                // Skip alignments which don't have at least min_similarity matches relative to larger node
                return;
            }

            // Update the graph
            alignment_graph.try_insert_node(id_a);
            alignment_graph.try_insert_node(id_b);

            alignment_graph.set_node_coverage(id_a, 0);
            alignment_graph.set_node_coverage(id_b, 0);

            alignment_graph.try_insert_edge(id_a, id_b, total_matches);

            alignment_graph.set_node_length(id_a, length_a);
            alignment_graph.set_node_length(id_b, length_b);

            cerr << "adding alignment: " << ordered_pair.first << ',' << ordered_pair.second << ',' << length_a << ',' << length_b << ',' << total_matches << '\n';
            cerr << "from graph: " << alignment_graph.get_node_length(id_a) << ',' << alignment_graph.get_node_length(id_b) << '\n';
        }
    });

    ContactGraph symmetrical_alignment_graph;

    bool symmetrical_edges_found = true;
    while (symmetrical_edges_found){
        sparse_hash_set <pair <int32_t, int32_t> > to_be_deleted;

        alignment_graph.for_each_edge([&](const pair<int32_t,int32_t> edge, int32_t weight){
            int32_t a = edge.first;
            int32_t b = edge.second;

            string a_name = id_map.get_name(a);
            string b_name = id_map.get_name(b);

            cerr << "Testing: " << a_name << ',' << b_name << '\n';

            int32_t a_best_neighbor = -1;
            int32_t a_best_value = -1;
            int32_t b_best_neighbor = -1;
            int32_t b_best_value = -1;

            cerr << '\t' << "-- a --" <<'\n';
            alignment_graph.for_each_node_neighbor(a, [&](int32_t other, const Node& n){
                auto w = alignment_graph.get_edge_weight(a, other);

                string other_name = id_map.get_name(other);
                cerr << '\t' << a_name << ',' << other_name << ',' << w << '\n';

                if (w > a_best_value){
                    cerr << "\tbest!" << '\n';
                    a_best_value = w;
                    a_best_neighbor = other;
                }
            });

            cerr << '\t' << "-- b --" <<'\n';
            alignment_graph.for_each_node_neighbor(b, [&](int32_t other, const Node& n){
                auto w = alignment_graph.get_edge_weight(b, other);

                string other_name = id_map.get_name(other);
                cerr << '\t' << b_name << ',' << other_name << ',' << w << '\n';

                if (w > b_best_value){
                    cerr << "\tbest!" << '\n';
                    b_best_value = w;
                    b_best_neighbor = other;
                }
            });

            // Cheap way to check existing overlap
            if (b_best_neighbor == a and a_best_neighbor == b){
                auto a_coverage = alignment_graph.get_node_coverage(a);
                auto b_coverage = alignment_graph.get_node_coverage(b);
                auto a_length = alignment_graph.get_node_length(a);
                auto b_length = alignment_graph.get_node_length(b);

                bool a_max = double(a_coverage) + double(weight) < double(a_length);
                bool b_max = double(b_coverage) + double(weight) < double(b_length);
                bool a_min = double(weight) > double(a_length)*0.2;
                bool b_min = double(weight) > double(b_length)*0.2;

                cerr << "current coverage on node a: " << a_coverage << " (length = " << a_length << ", weight = " << weight << ")" << '\n';
                cerr << "current coverage on node b: " << b_coverage << " (length = " << b_length << ", weight = " << weight << ")" << '\n';
                cerr << int(a_max) << ',' << int(b_max) << ',' << int(a_min) << ',' << int(b_min) << '\n';

                if (a_max and b_max and a_min and b_min){
                    symmetrical_alignment_graph.try_insert_node(a);
                    symmetrical_alignment_graph.try_insert_node(b);

                    symmetrical_alignment_graph.set_node_length(a, a_length);
                    symmetrical_alignment_graph.set_node_length(b, b_length);

                    symmetrical_alignment_graph.try_insert_edge(a, b, weight);

                    alignment_graph.increment_coverage(a,weight);
                    alignment_graph.increment_coverage(b,weight);

                    cerr << "good edge: " << a_name << ',' << b_name << '\n';
                    to_be_deleted.emplace(a,b);
                }
            }
        });

        symmetrical_edges_found = false;
        for (auto& [a,b]: to_be_deleted){
            alignment_graph.remove_edge(a,b);
            symmetrical_edges_found = true;
            cerr << "removing: " << id_map.get_name(a) << ',' << id_map.get_name(b) << '\n';
        }

        to_be_deleted.clear();
    }

    ofstream alignment_file(output_dir / "alignments.csv");
    alignment_file << "name_a" << ',' << "name_b" << ',' << "total_matches" << ',' << "symmetrical" << ',' << "color" << '\n';
    alignment_graph.for_each_edge([&](const pair<int32_t,int32_t> edge, int32_t weight){
        alignment_file << id_map.get_name(edge.first) << ',' << id_map.get_name(edge.second) << ',' << weight << ',' << 0 << ',' << "gray" << '\n';
        alignment_file << id_map.get_name(edge.second) << ',' << id_map.get_name(edge.first) << ',' << weight << ',' << 0 << ',' << "gray" << '\n';
    });

    symmetrical_alignment_graph.for_each_edge([&](const pair<int32_t,int32_t> edge, int32_t weight){
        auto a_length = alignment_graph.get_node_length(edge.first);
        auto b_length = alignment_graph.get_node_length(edge.second);

        pair <string, string> colors;

        if (a_length > b_length){
            colors = {"Cornflower Blue", "Tomato"};
        }
        else{
            colors = {"Tomato", "Cornflower Blue"};
        }

        alignment_file << id_map.get_name(edge.first) << ',' << id_map.get_name(edge.second) << ',' << weight << ',' << 1 << ',' << colors.first << '\n';
        alignment_file << id_map.get_name(edge.second) << ',' << id_map.get_name(edge.first) << ',' << weight << ',' << 1 << ',' << colors.second << '\n';
    });

    map<string,string> overlaps;

    hasher.get_symmetrical_matches(overlaps, 0.7);

    path output_path = output_dir / "pairs.csv";
    ofstream file(output_path);

    if (not file.good() or not file.is_open()){
        throw runtime_error("ERROR: could not write file: " + output_path.string());
    }

    file << "Name" << ',' << "Match" << ',' << "Color" << '\n';
    for (auto& [a,b]: overlaps){
        file << a << ',' << b << ',' << "Cornflower Blue" << '\n';
        file << b << ',' << a << ',' << "Tomato" << '\n';
    }

    cerr << t << "Done" << '\n';
}


int main (int argc, char* argv[]){
    path sam_path;
    path gfa_path;
    path output_dir;
    string required_prefix;
    int8_t min_mapq = 0;
    size_t n_threads = 1;

    CLI::App app{"App description"};

    app.add_option(
            "-g,--gfa",
            gfa_path,
            "Path to GFA containing assembly graph to be phased");

    app.add_option(
            "-o,--output_dir",
            output_dir,
            "Path to (nonexistent) directory where output will be stored")
            ->required();

    app.add_option(
            "-t,--threads",
            n_threads,
            "Maximum number of threads to use");

    CLI11_PARSE(app, argc, argv);

    phase_hic(output_dir, gfa_path, n_threads);

    return 0;
}


