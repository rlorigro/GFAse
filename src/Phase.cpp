#include "FixedBinarySequence.hpp"
#include "Bipartition.hpp"
#include "Phase.hpp"

#include "bdsg/overlays/packed_subgraph_overlay.hpp"

using bdsg::PackedSubgraphOverlay;

namespace gfase{


void get_chain_lengths(
        Bipartition& partition,
        IncrementalIdMap<string>& id_map,
        unordered_map<string, string>& diploid_path_names,
        map<size_t, size_t> chained_diploid_sizes,
        size_t rounding_unit
        ){

    // Collect some stats about the chains
    partition.for_each_subgraph([&](const HandleGraph& subgraph, size_t subgraph_index, bool partition){
        if (partition == 1){
            return;
        }

        size_t l_sum = 0;

        subgraph.for_each_handle([&](const handle_t& h){
            auto id = subgraph.get_id(h);
            auto name = id_map.get_name(id);
            size_t l = subgraph.get_length(h);

            // Average the lengths of the sides of each diploid bubble
            if (diploid_path_names.find(name) != diploid_path_names.end() and subgraph.get_node_count() > 1){
                l /= 2;
            }

            l_sum += l;
        });

        size_t s = l_sum/rounding_unit;
        chained_diploid_sizes[s] += 1;
    });
}


void phase_haplotype_paths(path gfa_path, size_t k, path paternal_kmers, path maternal_kmers, char path_delimiter) {
    HashGraph graph;
    IncrementalIdMap<string> id_map;
    KmerSets <FixedBinarySequence <uint64_t, 2> > ks(paternal_kmers, maternal_kmers);

    gfa_to_handle_graph(graph, id_map, gfa_path);

    plot_graph(graph, "start_graph");

    cerr << "Identifying diploid paths..." << '\n';

    unordered_map<string, string> diploid_path_names;
    unordered_set<string> haploid_path_names;
    find_diploid_paths(graph, diploid_path_names, haploid_path_names);

    cerr << "Extending paths by 1..." << '\n';

    vector <pair<path_handle_t, handle_t> > to_be_prepended;
    vector <pair<path_handle_t, handle_t> > to_be_appended;
    extend_paths(graph, to_be_prepended, to_be_appended);

    cerr << "Iterating path kmers..." << '\n';

    cerr << "\tNumber of components in graph: " << graph.get_path_count() << '\n';

    // Iterate paths and for each node, collect kmers if node is only covered by one path
    for (auto& [path_name, other_path_name]: diploid_path_names) {
        auto p = graph.get_path_handle(path_name);
        cerr << ">" << path_name << '\n';

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

    cerr << "Un-extending paths..." << '\n';

    un_extend_paths(graph, to_be_prepended, to_be_appended);
    to_be_prepended.clear();
    to_be_appended.clear();

    // ks.normalize_kmer_counts();

    // Open file and print header
    ofstream component_matrix_outfile("component_matrix_outfile.csv");
    component_matrix_outfile << "component_name,hap_0_paternal_count,hap_0_maternal_count,hap_1_paternal_count,hap_1_maternal_count \n"; 
    
    ks.for_each_component_matrix([&](const string& name, const array <array <float,2>, 2> component){
            component_matrix_outfile
                << component[0][KmerSets<string>::paternal_index] << ','
                << component[0][KmerSets<string>::maternal_index] << ','
                << component[1][KmerSets<string>::paternal_index] << ','
                << component[1][KmerSets<string>::maternal_index] << '\n';
    });

    // Debug
    ks.print_component_parent_conf_matrix();

    vector<HashGraph> connected_components;
    vector <IncrementalIdMap<string> > connected_component_ids;

    cerr << "Finding connected components..." << '\n';

    split_connected_components(graph, id_map, connected_components, connected_component_ids, true);

    cerr << "Unzipping..." << '\n';

    size_t rounding_unit = 10000;
    map<size_t, size_t> raw_diploid_sizes;
    map<size_t, size_t> chained_diploid_sizes;

    for (size_t c=0; c<connected_components.size(); c++){
        unzip(connected_components[c], connected_component_ids[c]);

        unordered_set<nid_t> diploid_nodes;
        auto& cc_graph = connected_components[c];
        auto& cc_id_map = connected_component_ids[c];

        // Collect some stats about the diploid nodes' lengths
//        cc_graph.for_each_handle([&](const handle_t& h){
//            auto id = cc_graph.get_id(h);
//            auto name = cc_id_map.get_name(id);
//            size_t l = cc_graph.get_length(h);
//
//            // Average the lengths of the sides of each diploid bubble
//            if (diploid_path_names.find(name) != diploid_path_names.end()){
//                l /= 2;
//            }
//
//            size_t s = l/rounding_unit;
//            raw_diploid_sizes[s] += 1;
//        });

        // Generate criteria for diploid node BFS
        cc_graph.for_each_path_handle([&](const path_handle_t& p){
            auto path_name = cc_graph.get_path_name(p);

            if (diploid_path_names.find(path_name) == diploid_path_names.end()){
                return true;
            }

            // Node names for haplotypes should match the paths that they were created from
            auto id = cc_id_map.get_id(path_name);

            // Add diploid nodes to the set
            diploid_nodes.emplace(id);

            return true;
        });

        Bipartition ploidy_bipartition(cc_graph, cc_id_map, diploid_nodes);
        ploidy_bipartition.partition();

        cerr << "PLOIDY subgraphs: " << '\n';
        ploidy_bipartition.print();

        unordered_set<nid_t> chain_nodes;
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
                        auto adjacent_subgraph_index = ploidy_bipartition.get_subgraph_index(id);

                        ploidy_bipartition.follow_subgraph_edges(adjacent_subgraph_index, false, [&](const handle_t& h2){
                            auto id2 = ploidy_bipartition.get_id(h2);
                            auto adjacent_subgraph_index2 = ploidy_bipartition.get_subgraph_index(id2);
                            second_degree_neighbors.emplace(adjacent_subgraph_index2);
                        });
                    });

                    ploidy_bipartition.follow_subgraph_edges(subgraph_index, false, [&](const handle_t& h){
                        auto id = ploidy_bipartition.get_id(h);
                        auto adjacent_subgraph_index = ploidy_bipartition.get_subgraph_index(id);

                        ploidy_bipartition.follow_subgraph_edges(adjacent_subgraph_index, true, [&](const handle_t& h2){
                            auto id2 = ploidy_bipartition.get_id(h2);
                            auto adjacent_subgraph_index2 = ploidy_bipartition.get_subgraph_index(id2);
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

        Bipartition chain_bipartition(cc_graph, cc_id_map, chain_nodes);
        chain_bipartition.partition();

        cerr << "CHAIN subgraphs: " << '\n';
        chain_bipartition.print();

        string filename_prefix = "component_" + to_string(c) + "_";
        ofstream file(filename_prefix + ".gfa");
        handle_graph_to_gfa(connected_components[c], connected_component_ids[c], file);

        ofstream test_gfa_meta(filename_prefix + "ploidy_metagraph.gfa");
        handle_graph_to_gfa(ploidy_bipartition.metagraph, test_gfa_meta);

        ofstream test_gfa_chain(filename_prefix + "chain_metagraph.gfa");
        handle_graph_to_gfa(chain_bipartition.metagraph, test_gfa_chain);

        ofstream test_csv_meta_chain(filename_prefix + "chain_metagraph.csv");
        chain_bipartition.write_meta_graph_csv(test_csv_meta_chain);

        ofstream test_csv_meta_ploidy(filename_prefix + "ploidy_metagraph.csv");
        ploidy_bipartition.write_meta_graph_csv(test_csv_meta_ploidy);

        ofstream test_csv_parent(filename_prefix + "chain_parent_graph.csv");
        chain_bipartition.write_parent_graph_csv(test_csv_parent);

        ofstream test_csv_parent_ploidy(filename_prefix + "ploidy_parent_graph.csv");
        ploidy_bipartition.write_parent_graph_csv(test_csv_parent_ploidy);

        ofstream test_csv_parent_chain(filename_prefix + "chain_parent_graph.csv");
        chain_bipartition.write_parent_graph_csv(test_csv_parent_chain);

        chain_bipartition.for_each_subgraph([&](const HandleGraph& subgraph, size_t subgraph_index, bool partition){
            queue <set <handle_t> > next_nodes;
            next_nodes.emplace();

            chain_bipartition.for_each_boundary_node_in_subgraph(subgraph_index, true, [&](const handle_t& h){
                next_nodes.back().emplace(h);
            });

            array <array <double, 2>, 2> matrix;

            while(not next_nodes.empty()){
                auto& nodes = next_nodes.front();

                if (nodes.size() == 1){
                    // Verify not diploid by name
                    auto name = cc_id_map.get_name(subgraph.get_id(*nodes.begin()));

                    if (diploid_path_names.count(name) > 0){
                        throw runtime_error("ERROR: singleton in chain is diploid: " + name);
                    }
                }
                else if (nodes.size() == 2){
                    // Verify nodes are diploid counterparts to one another
                    auto& node_a = *nodes.begin();
                    auto& node_b = *nodes.begin()++;

                    auto name_a = cc_id_map.get_name(subgraph.get_id(node_a));
                    auto name_b = cc_id_map.get_name(subgraph.get_id(node_b));

                    auto result = diploid_path_names.find(name_a);

                    if (result == diploid_path_names.end()){
                        throw runtime_error("ERROR: non diploid node in bubble of bubble chain: " + name_a);
                    }

                    if (result->second != name_b){
                        throw runtime_error("ERROR: nodes in bubble are not labeled as diploid counterparts: " + name_a + "," + name_b);
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
                    auto forward_score = matrix[haplotype_a][0] + matrix[haplotype_b][1];
                    auto flipped_score = matrix[haplotype_b][0] + matrix[haplotype_a][1];

                    if (forward_score > flipped_score){

                    }
                }
                else{
                    throw runtime_error("ERROR: diploid chain does not have 1 or 2 nodes in single position in chain");
                }

                next_nodes.pop();
            }
        });
    }

    cerr << "raw_diploid_sizes:" << '\n';
    for (auto& [size, count]: raw_diploid_sizes){
        cerr << size*rounding_unit << ',' << count << '\n';
    }

    cerr << "chained_diploid_sizes:" << '\n';
    for (auto& [size, count]: chained_diploid_sizes){
        cerr << size*rounding_unit << ',' << count << '\n';
    }


}

}