#include "IncrementalIdMap.hpp"
#include "gfa_to_handle.hpp"
#include "graph_utility.hpp"
#include "MultiContactGraph.hpp"
#include "Bipartition.hpp"
#include "hash_graph.hpp"
#include "Filesystem.hpp"
#include "Sequence.hpp"
#include "Hasher2.hpp"
#include "Timer.hpp"
#include "align.hpp"
#include "CLI11.hpp"
#include "Sam.hpp"
#include "Bam.hpp"
#include "minimap.h"

#include "SvgPlot.hpp"

using gfase::for_element_in_sam_file;
using gfase::random_multicontact_phase_search;
using gfase::unpaired_mappings_t;
using gfase::paired_mappings_t;
using gfase::contact_map_t;
using gfase::AlignmentBlock;
using gfase::AlignmentChain;

using gfase::gfa_to_handle_graph;
using gfase::IncrementalIdMap;
using gfase::MultiContactGraph;
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
#include <thread>

using std::unordered_map;
using std::thread;
using std::cref;
using std::ref;


void write_contact_map(
        path output_path,
        const MultiContactGraph& contact_graph,
        const IncrementalIdMap<string>& id_map){
    ofstream output_file(output_path);

    if (not output_file.is_open() or not output_file.good()){
        throw std::runtime_error("ERROR: could not write to file: " + output_path.string());
    }

    contact_graph.for_each_edge([&](const pair<int32_t,int32_t> edge, int32_t weight){
        output_file << id_map.get_name(edge.first) << ',' << id_map.get_name(edge.second) << ',' << weight << '\n';
    });
}


void write_graph_data(
        path output_directory,
        const Hasher2& hasher,
        const ContactGraph& contact_graph,
        const IncrementalIdMap<string>& id_map){

    path edges_path = output_directory / "edges.csv";
    path nodes_path = output_directory / "nodes.csv";

    ofstream edges_file(edges_path);
    ofstream nodes_file(nodes_path);

    if (not edges_file.is_open() or not edges_file.good()){
        throw std::runtime_error("ERROR: could not write to file: " + edges_path.string());
    }

    if (not nodes_file.is_open() or not nodes_file.good()){
        throw std::runtime_error("ERROR: could not write to file: " + nodes_path.string());
    }

    ContactGraph hash_graph;

    hasher.convert_to_contact_graph(hash_graph, id_map, 0.1, 10, 10);

    IncrementalIdMap<string> re_map(true);

    nodes_file << "id" << ',' << "name" << ',' << "length" << ',' << "contact_coverage" << ',' << "hash_coverage" << '\n';

    contact_graph.for_each_node([&](int32_t id, const Node& n){
        auto name = id_map.get_name(id);

        int64_t hash_coverage = 0;
        if (hash_graph.has_node(id)){
            hash_coverage = hash_graph.get_node_coverage(id);
        }

        auto new_id = re_map.try_insert(name);

        nodes_file << new_id << ',' << name << ',' << n.length << ',' << n.coverage << ',' << hash_coverage << '\n';
    });

    edges_file << "id_a" << ',' << "id_b" << ',' << "name_a" << ',' << "name_b" << ',' << "contact_weight" << ',' << "hash_weight" << '\n';

    // Iterate all contact graph edges and cross-reference with hash graph edges
    contact_graph.for_each_edge([&](const pair<int32_t,int32_t> edge, int32_t weight){
        auto& id_a = edge.first;
        auto& id_b = edge.second;

        // Edges will be ordered by their IDs {min(a,b), max(a,b)}
        auto name_a = id_map.get_name(id_a);
        auto name_b = id_map.get_name(id_b);

        auto new_id_a = re_map.get_id(name_a);
        auto new_id_b = re_map.get_id(name_b);

        auto hash_weight = hash_graph.get_edge_weight(id_a, id_b);
        auto contact_weight = contact_graph.get_edge_weight(id_a, id_b);

        edges_file << new_id_a << ',' << new_id_b << ',' << name_a << ',' << name_b << ',' << contact_weight << ',' << hash_weight << '\n';
    });
}


void update_contact_map(
        vector<SamElement>& alignments,
        MultiContactGraph& contact_graph,
        IncrementalIdMap<string>& id_map){

    // Iterate one triangle of the all-by-all matrix, adding up mapqs for reads on both end of the pair
    for (size_t i=0; i<alignments.size(); i++){
        auto& a = alignments[i];
        auto ref_id_a = int32_t(id_map.try_insert(a.ref_name));
        contact_graph.try_insert_node(ref_id_a, 0);

        contact_graph.increment_coverage(ref_id_a, 1);

        for (size_t j=i+1; j<alignments.size(); j++) {
            auto& b = alignments[j];
            auto ref_id_b = int32_t(id_map.try_insert(b.ref_name));
            contact_graph.try_insert_node(ref_id_b, 0);
            contact_graph.try_insert_edge(ref_id_a, ref_id_b);
            contact_graph.increment_edge_weight(ref_id_a, ref_id_b, 1);
        }
    }
}


void parse_unpaired_bam_file(
        path bam_path,
        MultiContactGraph& contact_graph,
        IncrementalIdMap<string>& id_map,
        string required_prefix,
        int8_t min_mapq){

    Bam reader(bam_path);

    size_t l = 0;
    string prev_query_name = "";
    vector<SamElement> alignments;

    reader.for_alignment_in_bam(false, [&](const SamElement& a){
        if (l == 0){
            prev_query_name = a.query_name;
        }

        if (prev_query_name != a.query_name){
            update_contact_map(alignments, contact_graph, id_map);
            alignments.clear();
        }

        // No information about reference contig, this alignment is unusable
        if (a.ref_name.empty()){
            return;
        }

        // Optionally filter by the contig names. E.g. "PR" in shasta
        bool valid_prefix = true;
        if (not required_prefix.empty()){
            for (size_t i=0; i<required_prefix.size(); i++){
                if (a.ref_name[i] != required_prefix[i]){
                    valid_prefix = false;
                    break;
                }
            }
        }

        if (valid_prefix) {
            // Only allow reads with mapq > min_mapq and not secondary
            if (a.mapq >= min_mapq and a.is_primary()) {
                alignments.emplace_back(a);
            }
        }

        l++;
        prev_query_name = a.query_name;
    });
}


void phase_contacts(
        const IncrementalIdMap<string>& id_map,
        MultiContactGraph& contact_graph,
        size_t n_threads
        ){

    vector<thread> threads;
    vector <pair <int32_t,int8_t> > best_partitions;
    vector<int32_t> ids;
    atomic<double> best_score = std::numeric_limits<double>::min();
    atomic<size_t> job_index = 0;
    mutex phase_mutex;
    size_t m_iterations = 10000;

    contact_graph.get_node_ids(ids);
    contact_graph.randomize_partitions();
    contact_graph.get_partitions(best_partitions);

    cerr << "start score: " << contact_graph.compute_total_consistency_score() << '\n';

    // Launch threads
    for (uint64_t i=0; i<n_threads; i++){
        try {
            threads.emplace_back(thread(
                    random_multicontact_phase_search,
                    contact_graph,
                    cref(ids),
                    ref(best_partitions),
                    ref(best_score),
                    ref(job_index),
                    ref(phase_mutex),
                    m_iterations
            ));
        } catch (const exception &e) {
            cerr << e.what() << "\n";
            exit(1);
        }
    }

    // Wait for threads to finish
    for (auto& t: threads){
        t.join();
    }

    contact_graph.set_partitions(best_partitions);
}


void phase_hic(path output_dir, path sam_path, path gfa_path, string required_prefix, int8_t min_mapq, size_t n_threads){
    Timer t;

    if (exists(output_dir)){
        throw runtime_error("ERROR: output directory exists already");
    }
    else {
        create_directories(output_dir);
    }

    // Check that input bam is readable before doing time-consuming steps
    {
        Bam reader(sam_path);
    }

    // Id-to-name bimap for reference contigs
    IncrementalIdMap<string> id_map(false);

    GfaReader reader(gfa_path);

    // TODO: move this into domain of bdsg graph instead of GFA reader
    string dummy_name;
    string dummy_seq;
    vector<Sequence> sequences = {Sequence(dummy_name,dummy_seq)};
    reader.for_each_sequence([&](string& name, string& sequence){
        id_map.insert(name);
        sequences.emplace_back(name, sequence);
    });

    // Hashing params
    double sample_rate = 0.04;
    size_t k = 22;
    size_t n_iterations = 6;

    // Only align the top n hits
    size_t max_hits = 5;

    // Hash results must have at least this percent similarity (A U B)/A, where A is larger.
    // Sequence lengths must be at least this ratio.
    // Resulting alignment coverage must be at least this amount on larger node.
    double min_similarity = 0.05;

    vector <pair <string,string> > to_be_aligned;

    get_alignment_candidates(
            sequences,
            id_map,
            to_be_aligned,
            output_dir,
            n_threads,
            sample_rate,
            k,
            n_iterations,
            max_hits,
            min_similarity
    );


    // TODO: convert edges to alts instead of using two parallel graphs... ?
    ContactGraph alignment_graph;
    ContactGraph symmetrical_alignment_graph;

    // Thread-related variables
    atomic<size_t> job_index = 0;
    vector<thread> threads;
    mutex output_mutex;

    // Launch threads
    for (uint64_t n=0; n<n_threads; n++){
        try {
            threads.emplace_back(thread(
                    construct_alignment_graph,
                    ref(to_be_aligned),
                    ref(sequences),
                    ref(id_map),
                    ref(alignment_graph),
                    min_similarity,
                    ref(output_mutex),
                    ref(job_index)
            ));
        } catch (const exception &e) {
            cerr << e.what() << "\n";
            exit(1);
        }
    }

    // Wait for threads to finish
    for (auto& n: threads){
        n.join();
    }

    get_best_overlaps(min_similarity, id_map, alignment_graph, symmetrical_alignment_graph);
    write_alignment_results_to_file(id_map, alignment_graph, symmetrical_alignment_graph, output_dir);

    cerr << t << "Done" << '\n';

    // Datastructures to represent linkages from hiC
    MultiContactGraph contact_graph;
    vector <vector <int32_t> > adjacency;

    cerr << t << "Loading alignments as contact map..." << '\n';

    if (sam_path.extension() == ".bam"){
        parse_unpaired_bam_file(sam_path, contact_graph, id_map, required_prefix, min_mapq);
    }
    else{
        throw runtime_error("ERROR: unrecognized extension for BAM input file: " + sam_path.extension().string());
    }

    for (auto& s: sequences){
        if (s.name.empty()){
            continue;
        }

        auto id = id_map.get_id(s.name);

        if (contact_graph.has_node(int32_t(id))) {
            contact_graph.set_node_length(int32_t(id), int32_t(s.size()));
        }
    }

    // Remove nodes that don't have any involvement in bubbles
    vector<int32_t> to_be_deleted;
    contact_graph.for_each_node([&](int32_t id){
        if (not symmetrical_alignment_graph.has_node(id)){
            to_be_deleted.emplace_back(id);
        }
    });

    for (auto& id: to_be_deleted){
        contact_graph.remove_node(id);
    }

    // Reallocate sparse hash map because something is going wrong during copying
    contact_graph.resize();

    // Add alts to graph
    symmetrical_alignment_graph.for_each_edge([&](const pair<int32_t,int32_t> edge, int32_t weight){
        auto [a,b] = edge;

        if (contact_graph.has_edge(a,b)){
            contact_graph.add_alt(a,b);
        }
    });

    phase_contacts(id_map, contact_graph, n_threads);

    cerr << t << "Writing phasing results to file... " << '\n';

    path contacts_output_path = output_dir / "contacts.csv";
    path phases_output_path = output_dir / "phases.csv";
    path config_output_path = output_dir / "config.csv";

    contact_graph.write_bandage_csv(phases_output_path, id_map);

    write_contact_map(contacts_output_path, contact_graph, id_map);

    // In lieu of actual chaining, dump some fasta files with the original GFA segments
    path phase_0_fasta_path = output_dir / "phase_0.fasta";
    ofstream phase_0_fasta(phase_0_fasta_path);

    path phase_1_fasta_path = output_dir / "phase_1.fasta";
    ofstream phase_1_fasta(phase_1_fasta_path);

    path unphased_fasta_path = output_dir / "unphased.fasta";
    ofstream unphased_fasta(unphased_fasta_path);

    if (not (phase_0_fasta.is_open() and phase_0_fasta.good())){
        throw runtime_error("ERROR: file could not be written: " + phase_0_fasta_path.string());
    }

    if (not (phase_1_fasta.is_open() and phase_1_fasta.good())){
        throw runtime_error("ERROR: file could not be written: " + phase_1_fasta_path.string());
    }

    if (not (unphased_fasta.is_open() and unphased_fasta.good())){
        throw runtime_error("ERROR: file could not be written: " + unphased_fasta_path.string());
    }

    // Iterate 1-based IDs and dump into fasta
    for (int64_t id=1; id<id_map.size()+1; id++){
        const auto& name = id_map.get_name(id);

        int8_t partition = 0;
        if (contact_graph.has_node(int32_t(id))){
            partition = contact_graph.get_partition(int32_t(id));
        }

        if (partition == 0){
            unphased_fasta << '>' << name << '\n';
            unphased_fasta << sequences.at(id).sequence << '\n';
        }
        else if (partition == -1){
            phase_0_fasta << '>' << name << '\n';
            phase_0_fasta << sequences.at(id).sequence << '\n';
        }
        else if (partition == 1){
            phase_1_fasta << '>' << name << '\n';
            phase_1_fasta << sequences.at(id).sequence << '\n';
        }
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
            "-i,--input",
            sam_path,
            "Path to SAM or BAM containing filtered, paired HiC reads")
            ->required();

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
            "-p,--prefix",
            required_prefix,
            "Prefix required in ref name for mapping to be counted");

    app.add_option(
            "-m,--min_mapq",
            min_mapq,
            "Minimum required mapq value for mapping to be counted");

    app.add_option(
            "-t,--threads",
            n_threads,
            "Maximum number of threads to use");

    CLI11_PARSE(app, argc, argv);

    phase_hic(output_dir, sam_path, gfa_path, required_prefix, min_mapq, n_threads);

    return 0;
}

