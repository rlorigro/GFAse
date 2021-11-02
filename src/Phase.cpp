#include "FixedBinarySequence.hpp"
#include "Bipartition.hpp"
#include "Phase.hpp"

#include "bdsg/overlays/packed_subgraph_overlay.hpp"

using bdsg::PackedSubgraphOverlay;

namespace gfase{


void count_kmers(
        const PathHandleGraph& graph,
        const IncrementalIdMap<string>& id_map,
        const unordered_map<string,string>& diploid_path_names,
        KmerSets <FixedBinarySequence <uint64_t,2> >& ks,
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
                FixedBinarySequence<uint64_t, 2> s(sequence);

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


void generate_ploidy_critera(
        const PathHandleGraph& graph,
        const IncrementalIdMap<string>& id_map,
        const unordered_map<string,string>& diploid_path_names,
        unordered_set<nid_t>& diploid_nodes
        ){

    graph.for_each_handle([&](const handle_t& h){
        auto name = id_map.get_name(graph.get_id(h));

        if (diploid_path_names.find(name) == diploid_path_names.end()){
            return true;
        }

        // Node names for haplotypes should match the paths that they were created from
        auto id = id_map.get_id(name);

        // Add diploid nodes to the set
        diploid_nodes.emplace(id);

        return true;
    });

}


void generate_chain_critera(
        Bipartition& ploidy_bipartition,
        unordered_set<nid_t>& chain_nodes
        ){

    ploidy_bipartition.for_each_subgraph([&](const HandleGraph& subgraph, size_t subgraph_index, bool partition){
        // Find singletons
        if (subgraph.get_node_count() == 1){
            if (partition == 0){
                subgraph.for_each_handle([&](const handle_t& h){
                    chain_nodes.emplace(subgraph.get_id(h));
                });
            }
            else{
                // If this is an unphased subgraph, check that it is not sharing its phased neighbors with any other
                // subgraphs, by doing a two-edge walk right/left and left/right
                unordered_set<size_t> second_degree_neighbors;

                ploidy_bipartition.follow_subgraph_edges(subgraph_index, true, [&](const handle_t& h){
                    auto id = ploidy_bipartition.get_id(h);
                    auto adjacent_subgraph_index = ploidy_bipartition.get_subgraph_index_of_parent_node(id);

                    ploidy_bipartition.follow_subgraph_edges(adjacent_subgraph_index, false, [&](const handle_t& h2){
                        auto id2 = ploidy_bipartition.get_id(h2);
                        auto adjacent_subgraph_index2 = ploidy_bipartition.get_subgraph_index_of_parent_node(id2);
                        second_degree_neighbors.emplace(adjacent_subgraph_index2);
                    });
                });

                ploidy_bipartition.follow_subgraph_edges(subgraph_index, false, [&](const handle_t& h){
                    auto id = ploidy_bipartition.get_id(h);
                    auto adjacent_subgraph_index = ploidy_bipartition.get_subgraph_index_of_parent_node(id);

                    ploidy_bipartition.follow_subgraph_edges(adjacent_subgraph_index, true, [&](const handle_t& h2){
                        auto id2 = ploidy_bipartition.get_id(h2);
                        auto adjacent_subgraph_index2 = ploidy_bipartition.get_subgraph_index_of_parent_node(id2);
                        second_degree_neighbors.emplace(adjacent_subgraph_index2);
                    });
                });

                // If there are no second degree neighbors, this unphased subgraph passes
                if (second_degree_neighbors.size() == 1){
                    ploidy_bipartition.for_each_handle_in_subgraph(subgraph_index, [&](const handle_t& h) {
                        chain_nodes.emplace(subgraph.get_id(h));
                    });
                }
            }
        }
        else if(subgraph.get_node_count() == 0){
            throw runtime_error("ERROR: subgraph in metagraph contains no nodes: " + to_string(subgraph_index));
        }
    });
}


void merge_diploid_singletons(const unordered_map<string,string>& diploid_path_names, Bipartition& chain_bipartition){
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
            auto result = diploid_path_names.find(singleton_name);

            if (result == diploid_path_names.end()){
                throw runtime_error("ERROR: diploid singleton has no counterpart in diploid names");
            }

            auto& other_name = result->second;
            auto other_id = chain_bipartition.get_id_of_parent_handle(other_name);

            auto other_subgraph_index = chain_bipartition.get_subgraph_index_of_parent_node(other_id);

            to_be_merged.emplace(min(subgraph_index,other_subgraph_index), max(subgraph_index,other_subgraph_index));
        }
    });

    for (auto& item: to_be_merged){
        chain_bipartition.merge_subgraphs(item.first, item.second);
    }
}


void phase_chains(
        Bipartition& chain_bipartition,
        MutablePathDeletableHandleGraph& graph,
        IncrementalIdMap<string>& id_map,
        const unordered_map<string,string>& diploid_path_names,
        const KmerSets <FixedBinarySequence <uint64_t,2> >& ks,
        unordered_set<string>& paternal_node_names,
        unordered_set<string>& maternal_node_names,
        ofstream& provenance_csv_file,
        string provenance_path_prefix,
        char path_delimiter
){

    chain_bipartition.for_each_subgraph([&](const HandleGraph& subgraph, size_t subgraph_index, bool partition){
        cerr << subgraph_index << '\n';

        // Skip unphased regions for now
        if (partition == 1){
            return;
        }

        array <array <double, 2>, 2> matrix;

        queue <set <handle_t> > q;
        set<handle_t> next_nodes;

        // Chain might terminate in an edge to another subgraph or in a dead end ("tip") so both need to be searched for
        set<handle_t> left_tips;
        set<handle_t> left_edge_nodes;

        // Find tips
        subgraph.for_each_handle([&](const handle_t& h){
            // Check the parent graph for edges
            if (graph.get_degree(h,true) == 0){
                left_tips.emplace(h);
            }
        });

        // Find edges to other subgraphs
        chain_bipartition.for_each_boundary_node_in_subgraph(subgraph_index, true, [&](const handle_t& h){
            left_edge_nodes.emplace(h);
            cerr << "queuing start node: " << id_map.get_name(graph.get_id(h)) << (graph.get_is_reverse(h) ? '-' : '+') << '\n';
        });

        // Make sure there are not both tips and edges
        if (left_tips.empty() and not left_edge_nodes.empty()){
            next_nodes = left_edge_nodes;
        }
        else if (not left_tips.empty() and left_edge_nodes.empty()){
            next_nodes = left_tips;
        }
        else if (not left_tips.empty() and not left_edge_nodes.empty()){
            throw runtime_error("ERROR: chain has left tips and left edge nodes (edges to other subgraph): " + to_string(subgraph_index));
        }

        path_handle_t maternal_path;
        path_handle_t paternal_path;

        if (not next_nodes.empty() and subgraph.get_node_count() > 1){
            q.emplace(next_nodes);

            string path_prefix = to_string(subgraph_index);

            string maternal_path_name = path_prefix + ".m";
            string paternal_path_name = path_prefix + ".p";

            maternal_path = graph.create_path_handle(maternal_path_name);
            paternal_path = graph.create_path_handle(paternal_path_name);

            paternal_node_names.emplace(paternal_path_name);
            maternal_node_names.emplace(maternal_path_name);
        }

        // Iterate each bubble or bridge and update the queue with the next nodes
        while(not q.empty()){
            auto& nodes = q.front();

            if (nodes.size() == 1){
                auto node = *nodes.begin();

                // Verify not diploid by name
                auto name = id_map.get_name(subgraph.get_id(node));

                // Handle case where singletons need phase assignment
                if (diploid_path_names.count(name) > 0){
                    throw runtime_error("ERROR: haploid node in non-singleton chain is flagged as diploid: " + name);
                }
                else {
                    // If it's a normal singleton just append it to both haplotype paths
                    graph.append_step(maternal_path, node);
                    graph.append_step(paternal_path, node);
                }
            }
            else if (nodes.size() == 2){
                // Verify nodes are diploid counterparts to one another
                auto& node_a = *nodes.begin();
                auto& node_b = *(++nodes.begin());

                auto name_a = id_map.get_name(subgraph.get_id(node_a));
                auto name_b = id_map.get_name(subgraph.get_id(node_b));

                auto result = diploid_path_names.find(name_a);

                if (result == diploid_path_names.end()){
                    throw runtime_error("ERROR: non diploid node in bubble of bubble chain: " + name_a);
                }

                if (result->second != name_b){
                    throw runtime_error("ERROR: nodes in bubble are not labeled as diploid counterparts: " + result->second + "," + name_b);
                }

                string component_name;
                size_t haplotype_a;
                size_t haplotype_b;

                // TODO: stop using names entirely!!
                // TODO: stop using names entirely!!
                // TODO: stop using names entirely!!
                tie(component_name, haplotype_a) = parse_path_string(name_a, path_delimiter);
                tie(component_name, haplotype_b) = parse_path_string(name_b, path_delimiter);

                // Fetch the kmer counts
                ks.get_matrix(component_name, matrix);

                // Forward score is the number of kmers supporting the orientation s.t. a == 0 and b == 1
                // Flipped score is the number of kmers supporting the orientation s.t. a == 1 and b == 0
                auto forward_score = matrix[haplotype_a][KmerSets<string>::paternal_index] + matrix[haplotype_b][KmerSets<string>::maternal_index];
                auto flipped_score = matrix[haplotype_b][KmerSets<string>::paternal_index] + matrix[haplotype_a][KmerSets<string>::maternal_index];

                // Choose the more supported orientation, defaulting to "forward orientation" if equal
                if (forward_score >= flipped_score){
                    graph.append_step(paternal_path, node_a);
                    graph.append_step(maternal_path, node_b);
                }
                else{
                    graph.append_step(maternal_path, node_a);
                    graph.append_step(paternal_path, node_b);
                }
            }
            else{
                cerr << "ERROR for subgraph_index: " << subgraph_index << " with nodes: " << '\n';
                for (auto& h: nodes){
                    cerr << '\t' << id_map.get_name(graph.get_id(h)) << '\n';
                }
                throw runtime_error("ERROR: diploid chain does not have 1 or 2 nodes in single position in chain");
            }

            // Find whatever comes next (if anything)
            next_nodes.clear();
            for (auto& node: nodes) {
                subgraph.follow_edges(node, false, [&](const handle_t& h){
                    next_nodes.emplace(h);
                });
            }

            if (not next_nodes.empty()){
                q.emplace(next_nodes);
            }

            q.pop();
        }
    });

    write_paths_to_csv(graph, id_map, provenance_csv_file, provenance_path_prefix);
    unzip(graph, id_map, false);
}


void phase_haplotype_paths(path gfa_path, size_t k, path paternal_kmers, path maternal_kmers, char path_delimiter) {
    HashGraph graph;
    IncrementalIdMap<string> id_map;
    KmerSets <FixedBinarySequence <uint64_t, 2> > ks(paternal_kmers, maternal_kmers);

    cerr << "Loading GFA..." << '\n';

    gfa_to_handle_graph(graph, id_map, gfa_path);

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
    component_matrix_outfile << "component_name,hap_0_paternal_count,hap_0_maternal_count,hap_1_paternal_count,hap_1_maternal_count \n"; 
    
    ks.for_each_component_matrix([&](const string& name, const array <array <double,2>, 2> matrix){
            component_matrix_outfile
                << matrix[0][KmerSets<string>::paternal_index] << ','
                << matrix[0][KmerSets<string>::maternal_index] << ','
                << matrix[1][KmerSets<string>::paternal_index] << ','
                << matrix[1][KmerSets<string>::maternal_index] << '\n';
    });

    // Debug
//    ks.print_component_parent_conf_matrix();

    vector<HashGraph> connected_components;
    vector <IncrementalIdMap<string> > connected_component_ids;

    cerr << "Finding connected components..." << '\n';

    split_connected_components(graph, id_map, connected_components, connected_component_ids, true);

    cerr << "Unzipping..." << '\n';

    ofstream maternal_fasta("maternal.fasta");
    ofstream paternal_fasta("paternal.fasta");
    ofstream unphased_fasta("unphased.fasta");

    path provenance_output_path = "phase_chains.csv";
    ofstream provenance_csv_file(provenance_output_path);

    if (not (provenance_csv_file.is_open() and provenance_csv_file.good())){
        throw runtime_error("ERROR: file could not be written: " + provenance_output_path.string());
    }

    provenance_csv_file << "path_name" << ',' << "n_steps" << ',' << "nodes" << '\n';

    for (size_t c=0; c<connected_components.size(); c++){
        unzip(connected_components[c], connected_component_ids[c], false);

        auto& cc_graph = connected_components[c];
        auto& cc_id_map = connected_component_ids[c];

        // Generate criteria for diploid node BFS
        unordered_set<nid_t> diploid_nodes;
        generate_ploidy_critera(cc_graph, cc_id_map, diploid_path_names, diploid_nodes);

        Bipartition ploidy_bipartition(cc_graph, cc_id_map, diploid_nodes);
        ploidy_bipartition.partition();

        // Generate criteria for node-chaining BFS
        unordered_set<nid_t> chain_nodes;
        generate_chain_critera(ploidy_bipartition, chain_nodes);

        Bipartition chain_bipartition(cc_graph, cc_id_map, chain_nodes);
        chain_bipartition.partition();

        string filename_prefix = "component_" + to_string(c) + "_";
        ofstream file(filename_prefix + ".gfa");
        handle_graph_to_gfa(connected_components[c], connected_component_ids[c], file);

        ofstream test_gfa_meta(filename_prefix + "ploidy_metagraph.gfa");
        handle_graph_to_gfa(ploidy_bipartition.metagraph, test_gfa_meta);

        ofstream test_csv_meta_ploidy(filename_prefix + "ploidy_metagraph.csv");
        ploidy_bipartition.write_meta_graph_csv(test_csv_meta_ploidy);

        ofstream test_csv_parent_ploidy(filename_prefix + "ploidy_parent_graph.csv");
        ploidy_bipartition.write_parent_graph_csv(test_csv_parent_ploidy);

        ofstream test_gfa_chain(filename_prefix + "chain_metagraph.gfa");
        handle_graph_to_gfa(chain_bipartition.metagraph, test_gfa_chain);

        ofstream test_csv_meta_chain(filename_prefix + "chain_metagraph.csv");
        chain_bipartition.write_meta_graph_csv(test_csv_meta_chain);

        ofstream test_csv_parent_chain(filename_prefix + "chain_parent_graph.csv");
        chain_bipartition.write_parent_graph_csv(test_csv_parent_chain);

        unordered_set<string> paternal_node_names;
        unordered_set<string> maternal_node_names;

        merge_diploid_singletons(diploid_path_names, chain_bipartition);

        ofstream test_gfa_chain_merged(filename_prefix + "chain_metagraph_merged.gfa");
        handle_graph_to_gfa(chain_bipartition.metagraph, test_gfa_chain_merged);

        ofstream test_csv_meta_chain_merged(filename_prefix + "chain_metagraph_merged.csv");
        chain_bipartition.write_meta_graph_csv(test_csv_meta_chain_merged);

        ofstream test_csv_parent_chain_merged(filename_prefix + "chain_parent_graph_merged.csv");
        chain_bipartition.write_parent_graph_csv(test_csv_parent_chain_merged);

        string component_path_prefix = to_string(c) + '.';
        phase_chains(
                chain_bipartition,
                cc_graph,
                cc_id_map,
                diploid_path_names,
                ks,
                paternal_node_names,
                maternal_node_names,
                provenance_csv_file,
                component_path_prefix,
                path_delimiter);

        ofstream test_gfa_phased(filename_prefix + "phased.gfa");
        handle_graph_to_gfa(cc_graph, cc_id_map, test_gfa_phased);

        cc_graph.for_each_handle([&](const handle_t& h){
            auto id = cc_graph.get_id(h);
            auto name = cc_id_map.get_name(id);

            if (maternal_node_names.count(name) > 0){
                maternal_fasta << '>' << component_path_prefix << name << '\n';
                maternal_fasta << cc_graph.get_sequence(h) << '\n';
            }
            else if (paternal_node_names.count(name) > 0){
                paternal_fasta << '>' << component_path_prefix << name << '\n';
                paternal_fasta << cc_graph.get_sequence(h) << '\n';
            }
            else {
                unphased_fasta << '>' << component_path_prefix << name << '\n';
                unphased_fasta << cc_graph.get_sequence(h) << '\n';
            }
        });
    }
}

}