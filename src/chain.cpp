#include "chain.hpp"

namespace gfase{

void generate_ploidy_criteria_from_bubble_graph(
        const PathHandleGraph& graph,
        const IncrementalIdMap<string>& id_map,
        const BubbleGraph& bubble_graph,
        unordered_set<nid_t>& diploid_nodes
){

    bubble_graph.for_each_node_id([&](const int32_t id){
        // Node names for haplotypes should match the paths that they were created from
        auto name = id_map.get_name(id);

        // Add diploid nodes to the set
        diploid_nodes.emplace(id);

        return true;
    });
}


void merge_diploid_singletons(const BubbleGraph& bubble_graph, Bipartition& chain_bipartition){
    unordered_set <pair <size_t,size_t> > to_be_merged;

    chain_bipartition.for_each_subgraph([&](const HandleGraph& subgraph, size_t subgraph_index, bool partition){
        // Look for phaseable subgraphs only (they contain at least one diploid node) and size == 1 (singleton)
        if (chain_bipartition.get_partition_of_subgraph(subgraph_index) == 0 and chain_bipartition.get_subgraph_size(subgraph_index) == 1){
            nid_t singleton_id;
            string singleton_name;

            chain_bipartition.for_each_handle_in_subgraph(subgraph_index, [&](const handle_t& h){
                singleton_id = chain_bipartition.get_id_of_parent_handle(h);
                singleton_name = chain_bipartition.get_name_of_parent_node(singleton_id);
            });

            // Find other diploid node and verify is also singleton
            nid_t other_id = bubble_graph.get_other_side(int32_t(singleton_id));

            size_t other_subgraph_index = -1;
            try {
                other_subgraph_index = chain_bipartition.get_subgraph_index_of_parent_node(other_id);
            }
            catch (exception& e){
                cerr << e.what() << '\n';
                cerr << "WARNING: skipping node because alt not found in graph: " << singleton_name << '\n';
                return;
            }

            // Use a defined ordering of singleton pairs to keep track of which have been visited
            to_be_merged.emplace(min(subgraph_index,other_subgraph_index), max(subgraph_index,other_subgraph_index));
        }
    });

    for (auto& item: to_be_merged){
        chain_bipartition.merge_subgraphs(item.first, item.second);
    }
}


void write_chaining_info_to_file(
        path output_dir,
        const Bipartition& ploidy_bipartition,
        const Bipartition& chain_bipartition,
        const PathHandleGraph& graph,
        const IncrementalIdMap<string>& id_map,
        const string& filename_prefix,
        size_t component_index,
        bool write_gfa
){

    size_t c = component_index;

    if (write_gfa) {
        path file_path = output_dir / "components" / to_string(c) / (filename_prefix + ".gfa");
        ofstream file(file_path);
        handle_graph_to_gfa(graph, id_map, file);
    }

    path test_gfa_chain_path = output_dir / "components" / to_string(c) / (filename_prefix + "chain_metagraph.gfa");
    ofstream test_gfa_chain(test_gfa_chain_path);

    if (not test_gfa_chain.is_open() or not test_gfa_chain.good()){
        throw runtime_error("ERROR: could not write to file: " + test_gfa_chain_path.string());
    }
    handle_graph_to_gfa(chain_bipartition.metagraph, test_gfa_chain);

    path test_csv_meta_chain_path = output_dir / "components" / to_string(c) / (filename_prefix + "chain_metagraph.csv");
    ofstream test_csv_meta_chain(test_csv_meta_chain_path);

    if (not test_csv_meta_chain.is_open() or not test_csv_meta_chain.good()){
        throw runtime_error("ERROR: could not write to file: " + test_csv_meta_chain_path.string());
    }
    chain_bipartition.write_meta_graph_csv(test_csv_meta_chain);

    path test_csv_parent_chain_path = output_dir / "components" / to_string(c) / (filename_prefix + "chain_parent_graph.csv");
    ofstream test_csv_parent_chain(test_csv_parent_chain_path);

    if (not test_csv_parent_chain.is_open() or not test_csv_parent_chain.good()){
        throw runtime_error("ERROR: could not write to file: " + test_csv_parent_chain_path.string());
    }
    chain_bipartition.write_parent_graph_csv(test_csv_parent_chain);
}


void chain_phased_gfa(
        MutablePathDeletableHandleGraph& graph,
        IncrementalIdMap<string>& id_map, const
        BubbleGraph& bubble_graph,
        path output_dir,
        bool write_gfa,
        bool write_fasta
        ){

    vector<HashGraph> connected_components;
    vector <IncrementalIdMap<string> > connected_component_ids;

    cerr << "Finding connected components..." << '\n';

    split_connected_components(graph, id_map, connected_components, true);

    cerr << "Chaining phased bubbles..." << '\n';

    for (size_t c=0; c<connected_components.size(); c++){
        unzip(connected_components[c], connected_component_ids[c], false);

        auto& cc_graph = connected_components[c];

        // Generate criteria for diploid node BFS
        unordered_set<nid_t> diploid_nodes;
        generate_ploidy_criteria_from_bubble_graph(cc_graph, id_map, bubble_graph, diploid_nodes);

        Bipartition ploidy_bipartition(cc_graph, id_map, diploid_nodes);
        ploidy_bipartition.partition();

        // Generate criteria for node-chaining BFS
        unordered_set<nid_t> chain_nodes;
        generate_chain_critera(ploidy_bipartition, chain_nodes);

        Bipartition chain_bipartition(cc_graph, id_map, chain_nodes);
        chain_bipartition.partition();

        create_directories(output_dir / "components" / to_string(c));

        string filename_prefix = "component_" + to_string(c) + "_";

        path test_gfa_meta_path = output_dir / "components" / to_string(c) / (filename_prefix + "ploidy_metagraph.gfa");
        ofstream test_gfa_meta(test_gfa_meta_path);

        if (not test_gfa_meta.is_open() or not test_gfa_meta.good()){
            throw runtime_error("ERROR: could not write to file: " + test_gfa_meta_path.string());
        }
        handle_graph_to_gfa(ploidy_bipartition.metagraph, test_gfa_meta);

        path test_csv_meta_ploidy_path = output_dir / "components" / to_string(c) / (filename_prefix + "ploidy_metagraph.csv");
        ofstream test_csv_meta_ploidy(test_csv_meta_ploidy_path);

        if (not test_csv_meta_ploidy.is_open() or not test_csv_meta_ploidy.good()){
            throw runtime_error("ERROR: could not write to file: " + test_csv_meta_ploidy_path.string());
        }
        ploidy_bipartition.write_meta_graph_csv(test_csv_meta_ploidy);

        path test_csv_parent_ploidy_path = output_dir / "components" / to_string(c) / (filename_prefix + "ploidy_parent_graph.csv");
        ofstream test_csv_parent_ploidy(test_csv_parent_ploidy_path);

        if (not test_csv_parent_ploidy.is_open() or not test_csv_parent_ploidy.good()){
            throw runtime_error("ERROR: could not write to file: " + test_csv_parent_ploidy_path.string());
        }
        ploidy_bipartition.write_parent_graph_csv(test_csv_parent_ploidy);

        write_chaining_info_to_file(output_dir, ploidy_bipartition, chain_bipartition, cc_graph, id_map, filename_prefix, c, write_gfa);

        merge_diploid_singletons(bubble_graph, chain_bipartition);

        path test_gfa_chain_merged_path = output_dir / "components" / to_string(c) / (filename_prefix + "chain_metagraph_merged.gfa");
        ofstream test_gfa_chain_merged(test_gfa_chain_merged_path);

        if (not test_gfa_chain_merged.is_open() or not test_gfa_chain_merged.good()){
            throw runtime_error("ERROR: could not write to file: " + test_gfa_chain_merged_path.string());
        }
        handle_graph_to_gfa(chain_bipartition.metagraph, test_gfa_chain_merged);

        path test_csv_meta_chain_merged_path = output_dir / "components" / to_string(c) / (filename_prefix + "chain_metagraph_merged.csv");
        ofstream test_csv_meta_chain_merged(test_csv_meta_chain_merged_path);

        if (not test_csv_meta_chain_merged.is_open() or not test_csv_meta_chain_merged.good()){
            throw runtime_error("ERROR: could not write to file: " + test_csv_meta_chain_merged_path.string());
        }
        chain_bipartition.write_meta_graph_csv(test_csv_meta_chain_merged);

        path test_csv_parent_chain_merged_path = output_dir / "components" / to_string(c) / (filename_prefix + "chain_parent_graph_merged.csv");
        ofstream test_csv_parent_chain_merged(test_csv_parent_chain_merged_path);

        if (not test_csv_parent_chain_merged.is_open() or not test_csv_parent_chain_merged.good()){
            throw runtime_error("ERROR: could not write to file: " + test_csv_parent_chain_merged_path.string());
        }
        chain_bipartition.write_parent_graph_csv(test_csv_parent_chain_merged);

        unordered_set<string> phase_0_node_names;
        unordered_set<string> phase_1_node_names;

        path_handle_t phase_0_path;
        path_handle_t phase_1_path;

        vector <vector <handle_t> > unphased_handles_per_component(connected_components.size());

        string component_path_prefix = to_string(c);
        auto prev_subgraph_index = numeric_limits<size_t>::max();
        for_element_in_bubble_chain(
                chain_bipartition,
                cc_graph,
                id_map,
                bubble_graph,
                [&](const vector<string>& node_names, size_t subgraph_index){

                    if (subgraph_index != prev_subgraph_index){
                        string path_prefix = component_path_prefix + '.' + to_string(subgraph_index);

                        string phase_0_path_name = path_prefix + ".0";
                        string phase_1_path_name = path_prefix + ".1";

                        phase_0_path = cc_graph.create_path_handle(phase_0_path_name);
                        phase_1_path = cc_graph.create_path_handle(phase_1_path_name);

                        phase_0_node_names.emplace(phase_0_path_name);
                        phase_1_node_names.emplace(phase_1_path_name);
                    }

                    // If there's only one node in the element, then it must be in-between bubbles
                    if (node_names.size() == 1){
                        auto node = cc_graph.get_handle(id_map.get_id(node_names[0]));

                        cc_graph.append_step(phase_0_path, node);
                        cc_graph.append_step(phase_1_path, node);
                    }
                    // If there are 2 nodes then it must be a bubble
                    else{
                        auto& name_a = node_names[0];
                        auto& name_b = node_names[1];

                        auto id_a = id_map.get_id(name_a);
                        auto id_b = id_map.get_id(name_b);

                        if (bubble_graph.find_bubble_id_of_node(int32_t(id_a)) != bubble_graph.find_bubble_id_of_node(int32_t(id_b))){
                            throw runtime_error("ERROR: nodes in diploid portion of chain are not of same bubble in bubble graph: " + (name_a + ',' + name_b));
                        }

                        // Only need one id to reach the corresponding bubble
                        auto bubble = bubble_graph.get_bubble_of_node(int32_t(id_a));

                        nid_t id_0 = bubble.first();
                        nid_t id_1 = bubble.second();

                        auto node_0 = cc_graph.get_handle(id_0);
                        auto node_1 = cc_graph.get_handle(id_1);

                        // Choose orientation supported by hi-c linkage (decided already, during phasing)
                        cc_graph.append_step(phase_0_path, node_0);
                        cc_graph.append_step(phase_1_path, node_1);
                    }

                    prev_subgraph_index = subgraph_index;
                });

        // Write out chaining paths before unzipping
        path provenance_csv_file_path = output_dir / "phase_chains.csv";
        ofstream provenance_csv_file(provenance_csv_file_path);

        if (not (provenance_csv_file.is_open() and provenance_csv_file.good())){
            throw runtime_error("ERROR: file could not be written: " + provenance_csv_file_path.string());
        }

        provenance_csv_file << "path_name" << ',' << "n_steps" << ',' << "nodes" << '\n';

        write_paths_to_csv(cc_graph, id_map, provenance_csv_file);

        unzip(cc_graph, id_map, false);

        if (write_gfa){
            // Write phased + unzipped gfa to file
            path test_gfa_phased_path = output_dir / "components" / to_string(c) / (filename_prefix + "phased.gfa");
            ofstream test_gfa_phased(test_gfa_phased_path);
            if (not test_gfa_phased.is_open() or not test_gfa_phased.good()){
                throw runtime_error("ERROR: could not write to file: " + test_gfa_phased_path.string());
            }

            handle_graph_to_gfa(cc_graph, id_map, test_gfa_phased);
        }


        if (write_fasta){
            path phase_0_fasta_path = output_dir / "phase_0.fasta";
            ofstream phase_0_fasta(phase_0_fasta_path);
            path phase_1_fasta_path = output_dir / "phase_1.fasta";
            ofstream phase_1_fasta(phase_1_fasta_path);
            path unphased_initial_fasta_path = output_dir / "unphased_initial.fasta";
            ofstream unphased_initial_fasta(unphased_initial_fasta_path);
            path unphased_fasta_path = output_dir / "unphased.fasta";
            ofstream unphased_fasta(unphased_fasta_path);

            if (not (phase_0_fasta.is_open() and phase_0_fasta.good())){
                throw runtime_error("ERROR: file could not be written: " + phase_0_fasta_path.string());
            }

            if (not (phase_1_fasta.is_open() and phase_1_fasta.good())){
                throw runtime_error("ERROR: file could not be written: " + phase_1_fasta_path.string());
            }

            if (not (unphased_initial_fasta.is_open() and unphased_initial_fasta.good())){
                throw runtime_error("ERROR: file could not be written: " + unphased_initial_fasta_path.string());
            }

            if (not (unphased_fasta.is_open() and unphased_fasta.good())){
                throw runtime_error("ERROR: file could not be written: " + unphased_fasta_path.string());
            }

            cc_graph.for_each_handle([&](const handle_t& h){
                auto id = cc_graph.get_id(h);
                auto name = id_map.get_name(id);

                if (phase_0_node_names.count(name) > 0){
                    phase_0_fasta << '>' << name << '\n';
                    phase_0_fasta << cc_graph.get_sequence(h) << '\n';
                }
                else if (phase_1_node_names.count(name) > 0){
                    phase_1_fasta << '>' << name << '\n';
                    phase_1_fasta << cc_graph.get_sequence(h) << '\n';
                }
                else {
                    unphased_initial_fasta << '>' << name << '\n';
                    unphased_initial_fasta << cc_graph.get_sequence(h) << '\n';
                    unphased_handles_per_component[c].emplace_back(h);
                }
            });

        }
    }
}

}