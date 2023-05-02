#ifndef GFASE_PHASE_HPP
#define GFASE_PHASE_HPP

#include "FixedBinarySequence.hpp"
#include "HaplotypePathKmer.hpp"
#include "IncrementalIdMap.hpp"
#include "gfa_to_handle.hpp"
#include "handle_to_gfa.hpp"
#include "graph_utility.hpp"
#include "Bipartition.hpp"
#include "BubbleGraph.hpp"
#include "Filesystem.hpp"
#include "GfaReader.hpp"
#include "KmerSets.hpp"
#include "Color.hpp"
#include "CLI11.hpp"

#include "bdsg/hash_graph.hpp"
#include "bdsg/overlays/packed_subgraph_overlay.hpp"

#include <string>

using ghc::filesystem::path;

using bdsg::HashGraph;
using bdsg::MutablePathMutableHandleGraph;
using bdsg::MutablePathDeletableHandleGraph;
using bdsg::PackedSubgraphOverlay;
using handlegraph::path_handle_t;
using handlegraph::step_handle_t;
using handlegraph::handle_t;

using std::string;
using std::cout;
using std::cerr;


namespace gfase {

//// Find any nodes that are adjacent to the beginning and end of a path, as long as they are the only adjacent node
//pair<handle_t, bool> find_singleton_adjacent_handle(const PathHandleGraph& graph, const handle_t& h, bool left);
//
//
//void extend_paths(MutablePathMutableHandleGraph& graph);


//template <class T, size_t T2> void phase(
//        path gfa_path,
//        size_t k,
//        path paternal_kmers,
//        path maternal_kmers,
//        char path_delimiter='.');
//
//
//void phase_k(path gfa_path, size_t k, path paternal_kmers, path maternal_kmers, char path_delimiter);


template <class T, size_t T2> void count_kmers(
        const PathHandleGraph& graph,
        const IncrementalIdMap<string>& id_map,
        const unordered_map<string,string>& diploid_path_names,
        KmerSets <FixedBinarySequence <T,T2> >& ks,
        size_t k,
        char path_delimiter) {

    // Iterate paths and for each node, collect kmers if node is only covered by one path
    for (auto& [path_name, other_path_name]: diploid_path_names) {
        auto p = graph.get_path_handle(path_name);
        HaplotypePathKmer kmer(graph, p, k);

        string component_name;
        size_t haplotype;

        // TODO: stop using names entirely!!
        // TODO: stop using names entirely!!
        // TODO: stop using names entirely!!
        tie(component_name, haplotype) = parse_path_string(path_name, path_delimiter);

        kmer.for_each_haploid_kmer([&](const deque<char>& sequence){
            try {
                FixedBinarySequence<T, T2> s(sequence);

                // Compare kmer to parental kmers
                ks.increment_parental_kmer_count(component_name, haplotype, s);
            }
            catch(exception& e){
                auto node_name = id_map.get_name(graph.get_id(graph.get_handle_of_step(kmer.get_step_of_kmer_end())));
                cerr << e.what() << '\n';
                throw runtime_error("Error parsing sequence for node: " + node_name);
            }
        });
    }
}


void generate_ploidy_critera_from_bubble_bimap(
        const PathHandleGraph& graph,
        const IncrementalIdMap<string>& id_map,
        const unordered_map<string,string>& diploid_path_names,
        unordered_set<nid_t>& diploid_nodes
);

void generate_ploidy_criteria_from_bubble_graph(
        const PathHandleGraph& graph,
        const IncrementalIdMap<string>& id_map,
        const unordered_map<string,string>& diploid_path_names,
        unordered_set<nid_t>& diploid_nodes
);

// TODO: convert this side of the pipeline (trio phasing) to use BubbleGraph instead of "diploid_path_names" bimap
void merge_diploid_singletons(const unordered_map<string,string>& diploid_path_names, Bipartition& chain_bipartition);


template <class T, size_t T2> void order_bubble(Bubble<string>& b, const KmerSets<FixedBinarySequence <T, T2> >& ks){
    array <array <double, 2>, 2> matrix;

    string component_name;
    size_t haplotype_a;
    size_t haplotype_b;

    // TODO: stop using names entirely!!
    // TODO: stop using names entirely!!
    // TODO: stop using names entirely!!
    tie(component_name, haplotype_a) = parse_path_string(b.first(), ks.path_delimiter);
    tie(component_name, haplotype_b) = parse_path_string(b.second(), ks.path_delimiter);

    // Fetch the kmer counts
    ks.get_matrix(component_name, matrix);

    // Forward score is the number of kmers supporting the orientation s.t. a == 0 and b == 1
    // Flipped score is the number of kmers supporting the orientation s.t. a == 1 and b == 0
    auto forward_score = matrix[haplotype_a][KmerSets<string>::paternal_index] + matrix[haplotype_b][KmerSets<string>::maternal_index];
    auto flipped_score = matrix[haplotype_b][KmerSets<string>::paternal_index] + matrix[haplotype_a][KmerSets<string>::maternal_index];

    if (flipped_score > forward_score){
        b.flip();
    }
}


// TODO: delete this method in favor of bubble graph method
void for_element_in_bubble_chain(
        Bipartition& chain_bipartition,
        MutablePathDeletableHandleGraph& graph,
        IncrementalIdMap<string>& id_map,
        const unordered_map<string,string>& diploid_path_names,
        const function<void(const vector<string>& node_names, size_t subgraph_index)>& f
);


template <class T, size_t T2>
void phase(path gfa_path, size_t k, path paternal_kmers, path maternal_kmers, path output_directory, char path_delimiter) {
    if (exists(output_directory)){
        throw runtime_error("ERROR: output directory exists already");
    }
    else {
        create_directories(output_directory);
        create_directories(output_directory / "components");
    }

    size_t min_unphased_contig_length = 100000;
    double min_total_kmers = 30;
    double score_threshold = 0.02;

    HashGraph graph;
    IncrementalIdMap<string> id_map;
    Overlaps overlaps;
    KmerSets <FixedBinarySequence <T, T2> > ks(paternal_kmers, maternal_kmers);

    if (ks.get_k() != k){
        throw runtime_error("ERROR: kmers in file " + to_string(ks.get_k()) + " do not match k " + to_string(k));
    }

    cerr << "Loading GFA..." << '\n';

    gfa_to_handle_graph(graph, id_map, overlaps, gfa_path, false);

    cerr << "\tNumber of components in graph: " << graph.get_path_count() << '\n';

    cerr << "Identifying diploid paths..." << '\n';

    unordered_map<string, string> diploid_path_names;
    unordered_set<string> haploid_path_names;
    find_diploid_paths(graph, diploid_path_names, haploid_path_names);

    cerr << "Extending paths by 1..." << '\n';

    vector <pair<path_handle_t, handle_t> > to_be_prepended;
    vector <pair<path_handle_t, handle_t> > to_be_appended;
    extend_paths(graph, to_be_prepended, to_be_appended);

    cerr << "Iterating path kmers..." << '\n';

    count_kmers(graph, id_map, diploid_path_names, ks, k, path_delimiter);

    cerr << "Un-extending paths..." << '\n';

    un_extend_paths(graph, to_be_prepended, to_be_appended);
    to_be_prepended.clear();
    to_be_appended.clear();

    // Open file and print header
    ofstream component_matrix_outfile("component_matrix_outfile.csv");
    component_matrix_outfile << "component_name,hap_0_paternal_count,hap_0_maternal_count,hap_1_paternal_count,hap_1_maternal_count\n";

    ks.for_each_component_matrix([&](const string& name, const array <array <double,2>, 2> matrix){
        component_matrix_outfile
                << name << ','
                << matrix[0][KmerSets<string>::paternal_index] << ','
                << matrix[0][KmerSets<string>::maternal_index] << ','
                << matrix[1][KmerSets<string>::paternal_index] << ','
                << matrix[1][KmerSets<string>::maternal_index] << '\n';
    });

    // Debug
//    ks.print_component_parent_conf_matrix();

    vector<HashGraph> connected_components;
    vector <IncrementalIdMap<string> > connected_component_ids;
    vector<Overlaps> connected_component_overlaps;

    cerr << "Finding connected components..." << '\n';

    split_connected_components(graph, id_map, overlaps, connected_components, connected_component_ids, connected_component_overlaps, true);

    cerr << "Unzipping..." << '\n';

    ofstream maternal_fasta(output_directory / "maternal.fasta");
    ofstream paternal_fasta(output_directory / "paternal.fasta");
    ofstream unphased_initial_fasta(output_directory / "unphased_initial.fasta");
    ofstream unphased_fasta(output_directory / "unphased.fasta");

    path provenance_output_path = output_directory / "phase_chains.csv";
    ofstream provenance_csv_file(provenance_output_path);

    if (not (provenance_csv_file.is_open() and provenance_csv_file.good())){
        throw runtime_error("ERROR: file could not be written: " + provenance_output_path.string());
    }

    provenance_csv_file << "path_name" << ',' << "n_steps" << ',' << "nodes" << '\n';

    vector <vector <handle_t> > unphased_handles_per_component(connected_components.size());

    for (size_t c=0; c<connected_components.size(); c++){
        
        auto& cc_graph = connected_components[c];
        auto& cc_id_map = connected_component_ids[c];
        auto& cc_overlaps = connected_component_overlaps[c];
        
        unzip(cc_graph, cc_id_map, cc_overlaps, false);

        // Generate criteria for diploid node BFS
        unordered_set<nid_t> diploid_nodes;
        generate_ploidy_critera_from_bubble_bimap(cc_graph, cc_id_map, diploid_path_names, diploid_nodes);

        Bipartition ploidy_bipartition(cc_graph, cc_id_map, diploid_nodes);
        ploidy_bipartition.partition();

        // Generate criteria for node-chaining BFS
        unordered_set<nid_t> chain_nodes;
        generate_chain_critera(ploidy_bipartition, chain_nodes);

        Bipartition chain_bipartition(cc_graph, cc_id_map, chain_nodes);
        chain_bipartition.partition();

        string filename_prefix = output_directory / "components" / ("component_" + to_string(c) + "_");
        ofstream file(filename_prefix + ".gfa");
        handle_graph_to_gfa(cc_graph, cc_id_map, cc_overlaps, file);

        ofstream test_gfa_meta(filename_prefix + "ploidy_metagraph.gfa");
        handle_graph_to_gfa(ploidy_bipartition.metagraph, test_gfa_meta, Overlaps());

        ofstream test_csv_meta_ploidy(filename_prefix + "ploidy_metagraph.csv");
        ploidy_bipartition.write_meta_graph_csv(test_csv_meta_ploidy);

        ofstream test_csv_parent_ploidy(filename_prefix + "ploidy_parent_graph.csv");
        ploidy_bipartition.write_parent_graph_csv(test_csv_parent_ploidy);

        ofstream test_gfa_chain(filename_prefix + "chain_metagraph.gfa");
        handle_graph_to_gfa(chain_bipartition.metagraph, test_gfa_chain, Overlaps());

        ofstream test_csv_meta_chain(filename_prefix + "chain_metagraph.csv");
        chain_bipartition.write_meta_graph_csv(test_csv_meta_chain);

        ofstream test_csv_parent_chain(filename_prefix + "chain_parent_graph.csv");
        chain_bipartition.write_parent_graph_csv(test_csv_parent_chain);

        unordered_set<string> paternal_node_names;
        unordered_set<string> maternal_node_names;

        merge_diploid_singletons(diploid_path_names, chain_bipartition);

        ofstream test_gfa_chain_merged(filename_prefix + "chain_metagraph_merged.gfa");
        handle_graph_to_gfa(chain_bipartition.metagraph, test_gfa_chain_merged, Overlaps());

        ofstream test_csv_meta_chain_merged(filename_prefix + "chain_metagraph_merged.csv");
        chain_bipartition.write_meta_graph_csv(test_csv_meta_chain_merged);

        ofstream test_csv_parent_chain_merged(filename_prefix + "chain_parent_graph_merged.csv");
        chain_bipartition.write_parent_graph_csv(test_csv_parent_chain_merged);

        string component_path_prefix = to_string(c);

        path_handle_t maternal_path;
        path_handle_t paternal_path;

        auto prev_subgraph_index = numeric_limits<size_t>::max();
        for_element_in_bubble_chain(
                chain_bipartition,
                cc_graph,
                cc_id_map,
                diploid_path_names,
                [&](const vector<string>& node_names, size_t subgraph_index){

            if (subgraph_index != prev_subgraph_index){
                string path_prefix = component_path_prefix + '.' + to_string(subgraph_index);

                string maternal_path_name = path_prefix + ".m";
                string paternal_path_name = path_prefix + ".p";

                maternal_path = cc_graph.create_path_handle(maternal_path_name);
                paternal_path = cc_graph.create_path_handle(paternal_path_name);

                paternal_node_names.emplace(paternal_path_name);
                maternal_node_names.emplace(maternal_path_name);
            }

            // If there's only one node in the element, then it must be in-between bubbles
            if (node_names.size() == 1){
                auto node = cc_graph.get_handle(cc_id_map.get_id(node_names[0]));

                cc_graph.append_step(maternal_path, node);
                cc_graph.append_step(paternal_path, node);
            }
            // If there are 2 nodes then it must be a bubble
            else{
                auto& name_a = node_names[0];
                auto& name_b = node_names[1];

                Bubble<string> b(name_a, name_b, 0);
                order_bubble(b, ks);

                auto node_a = cc_graph.get_handle(cc_id_map.get_id(name_a));
                auto node_b = cc_graph.get_handle(cc_id_map.get_id(name_b));

                // Choose the more supported orientation, defaulting to "forward orientation" if equal
                if (b.phase == 0){
                    cc_graph.append_step(paternal_path, node_a);
                    cc_graph.append_step(maternal_path, node_b);
                }
                else{
                    cc_graph.append_step(maternal_path, node_a);
                    cc_graph.append_step(paternal_path, node_b);
                }
            }

            prev_subgraph_index = subgraph_index;
        });

        unzip(cc_graph, cc_id_map, cc_overlaps, false);
        write_paths_to_csv(cc_graph, cc_id_map, provenance_csv_file);

        ofstream test_gfa_phased(output_directory / (filename_prefix + "phased.gfa"));
        handle_graph_to_gfa(cc_graph, cc_id_map, cc_overlaps, test_gfa_phased);

        cc_graph.for_each_handle([&](const handle_t& h){
            auto id = cc_graph.get_id(h);
            auto name = cc_id_map.get_name(id);

            if (maternal_node_names.count(name) > 0){
                maternal_fasta << '>' << name << '\n';
                maternal_fasta << cc_graph.get_sequence(h) << '\n';
            }
            else if (paternal_node_names.count(name) > 0){
                paternal_fasta << '>' << name << '\n';
                paternal_fasta << cc_graph.get_sequence(h) << '\n';
            }
            else {
                unphased_initial_fasta << '>' << name << '\n';
                unphased_initial_fasta << cc_graph.get_sequence(h) << '\n';
                unphased_handles_per_component[c].emplace_back(h);
            }
        });
    }

    ofstream unphased_parental_counts(output_directory / "unphased_parental_counts.csv");
    unphased_parental_counts << "name" << ',' << "maternal_count" << ',' << "paternal_count" << ',' << "unique_maternal_count" << ',' << "unique_paternal_count" << ',' << "score" << ',' << "unique_score" << ',' << "color" << '\n';

    ofstream unphased_kmers_log(output_directory / "unphased_kmers.csv");
    unphased_kmers_log << "count" << ',' << "is_mat" << ',' << "is_pat" << '\n';

    sparse_hash_map <FixedBinarySequence<T,T2>, size_t> unphased_kmer_counts;
    Coolwarm colormap;

    // Count k-mers, so we can later eliminate non-unique ones
    for (size_t c=0; c<connected_components.size(); c++){
        auto& cc_graph = connected_components[c];

        for (auto& h: unphased_handles_per_component[c]){
            size_t contig_length = cc_graph.get_length(h);

            string sequence = cc_graph.get_sequence(h);
            string initial_kmer = sequence.substr(0,k);
            FixedBinarySequence<T,T2> kmer(initial_kmer);

            unphased_kmer_counts[kmer]++;

            for (size_t i = k; i<contig_length; i++){
                kmer.shift(sequence[i], k);
                unphased_kmer_counts[kmer]++;
            }
        }
    }

    // Eliminate non-unique kmers
    sparse_hash_set <FixedBinarySequence<T,T2> > unique_kmers;
    for (auto& [kmer, count]: unphased_kmer_counts){
        if (count == 1){
            unique_kmers.emplace(kmer);
        }

        bool is_mat = ks.is_maternal(kmer);
        bool is_pat = ks.is_paternal(kmer);

        unphased_kmers_log << count << ',' << is_mat << ',' << is_pat << '\n';
    }

    unphased_kmer_counts.clear();

    // Handle unphased nodes
    for (size_t c=0; c<connected_components.size(); c++){
        auto& cc_graph = connected_components[c];
        auto& cc_id_map = connected_component_ids[c];

        for (auto& h: unphased_handles_per_component[c]){
            size_t contig_length = cc_graph.get_length(h);

            auto name = cc_id_map.get_name(cc_graph.get_id(h));

            string sequence = cc_graph.get_sequence(h);
            string initial_kmer = sequence.substr(0,k);
            FixedBinarySequence<T,T2> kmer(initial_kmer);

            double maternal_count = ks.is_maternal(kmer);
            double paternal_count = ks.is_paternal(kmer);

            double unique_maternal_count = ks.is_maternal(kmer) and unique_kmers.contains(kmer);
            double unique_paternal_count = ks.is_paternal(kmer) and unique_kmers.contains(kmer);

            for (size_t i = k; i<cc_graph.get_length(h); i++){
                kmer.shift(sequence[i], k);

                bool is_mat = ks.is_maternal(kmer);
                bool is_pat = ks.is_paternal(kmer);

                maternal_count += is_mat;
                paternal_count += is_pat;

                if (unique_kmers.contains(kmer)) {
                    unique_maternal_count += is_mat;
                    unique_paternal_count += is_pat;
                }
            }

            double normalized_paternal_score = (paternal_count + 0.000001)/ks.n_paternal_kmers();
            double normalized_maternal_score = (maternal_count + 0.000001)/ks.n_maternal_kmers();

            double score = log2(normalized_paternal_score/normalized_maternal_score);
            double unique_score = log2((unique_paternal_count + 0.000001)/(unique_maternal_count + 0.000001));

            double color_index = 0.5 + 0.5*(score/3);   // Saturate at 3x

            color_index = min(1.0,color_index);
            color_index = max(0.0,color_index);

            auto color = colormap.get_rgb(color_index);
            string hex_color = rgb_to_hex(color[0], color[1], color[2]);

            unphased_parental_counts << name << ',' << maternal_count << ',' << paternal_count << ',' << unique_maternal_count << ',' << unique_paternal_count << ',' << score << ',' << unique_score << ',' << '#' << hex_color << '\n';

            if ((contig_length > min_unphased_contig_length) and (maternal_count + paternal_count > min_total_kmers)){
                if (unique_score > score_threshold){
                    paternal_fasta << '>' << name << '\n';
                    paternal_fasta << cc_graph.get_sequence(h) << '\n';
                }
                else if (unique_score < -score_threshold){
                    maternal_fasta << '>' << name << '\n';
                    maternal_fasta << cc_graph.get_sequence(h) << '\n';
                }
                else{
                    unphased_fasta << '>' << name << '\n';
                    unphased_fasta << cc_graph.get_sequence(h) << '\n';
                }
            }
            else{
                // These contigs don't meet minimum criteria to even evaluate their score
                unphased_fasta << '>' << name << '\n';
                unphased_fasta << cc_graph.get_sequence(h) << '\n';
            }
        }
    }
}


void phase_k(path gfa_path, size_t k, path paternal_kmers, path maternal_kmers, path output_directory, char path_delimiter='.');


}

#endif //GFASE_PHASE_HPP
