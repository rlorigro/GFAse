#include "Phase.hpp"

namespace gfase{



void generate_ploidy_critera_from_bubble_bimap(
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


void generate_ploidy_criteria_from_bubble_graph(
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


// TODO: convert this side of the pipeline (trio phasing) to use BubbleGraph instead of "diploid_path_names" bimap
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

            // Use a defined ordering of singleton pairs to keep track of which have been visited
            to_be_merged.emplace(min(subgraph_index,other_subgraph_index), max(subgraph_index,other_subgraph_index));
        }
    });

    for (auto& item: to_be_merged){
        chain_bipartition.merge_subgraphs(item.first, item.second);
    }
}


// TODO: delete this method in favor of bubble graph method
void for_element_in_bubble_chain(
        Bipartition& chain_bipartition,
        MutablePathDeletableHandleGraph& graph,
        IncrementalIdMap<string>& id_map,
        const unordered_map<string,string>& diploid_path_names,
        const function<void(const vector<string>& node_names, size_t subgraph_index)>& f
){

    chain_bipartition.for_each_subgraph([&](const HandleGraph& subgraph, size_t subgraph_index, bool partition){
        cerr << subgraph_index << '\n';

        // Skip unphased regions for now
        if (partition == 1){
            return;
        }

        queue <set <handle_t> > q;
        set<nid_t> visited;
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

//        path_handle_t maternal_path;
//        path_handle_t paternal_path;

        // Initialize things for this chain
        if (not next_nodes.empty() and subgraph.get_node_count() > 1){
            q.emplace(next_nodes);

//            string path_prefix = provenance_path_prefix + '.' + to_string(subgraph_index);
//
//            string maternal_path_name = path_prefix + ".m";
//            string paternal_path_name = path_prefix + ".p";
//
//            maternal_path = graph.create_path_handle(maternal_path_name);
//            paternal_path = graph.create_path_handle(paternal_path_name);
//
//            phase_0_node_names.emplace(paternal_path_name);
//            phase_1_node_names.emplace(maternal_path_name);
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
                    vector<string> n = {name};
                    f(n, subgraph_index);
//                    // If it's a normal singleton just append it to both haplotype paths
//                    graph.append_step(maternal_path, node);
//                    graph.append_step(paternal_path, node);
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

                vector<string> n = {name_a, name_b};
                f(n, subgraph_index);

//                Bubble<string> b(name_a, name_b, 0);
//                order_bubble(b);
//
//                // Choose the more supported orientation, defaulting to "forward orientation" if equal
//                if (b.phase == 0){
//                    graph.append_step(paternal_path, node_a);
//                    graph.append_step(maternal_path, node_b);
//                }
//                else{
//                    graph.append_step(maternal_path, node_a);
//                    graph.append_step(paternal_path, node_b);
//                }
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
                    auto h_id = subgraph.get_id(h);

                    // Avoid adding self-looped or reversing nodes more than once to the queue
                    if (visited.count(h_id) == 0) {
                        next_nodes.emplace(h);
                        visited.emplace(h_id);
                    }
                });
            }

            if (not next_nodes.empty()){
                q.emplace(next_nodes);
            }

            q.pop();
        }
    });
}


void phase_k(path gfa_path, size_t k, path paternal_kmers, path maternal_kmers, char path_delimiter){
    if (k < 4){
        throw runtime_error("ERROR: must choose a k value larger than 4");
    }
        // Min = 8 bits, max = 16 bits
    else if (k >= 4 and k <= 8){
        phase<uint16_t,1>(gfa_path, k, paternal_kmers, maternal_kmers, path_delimiter);
    }
        // Min = 18 bits, max = 24 bits
    else if (k > 8 and k <= 12){
        phase<uint8_t,3>(gfa_path, k, paternal_kmers, maternal_kmers, path_delimiter);
    }
        // Min = 26 bits, max = 32 bits
    else if (k > 12 and k <= 16){
        phase<uint32_t,1>(gfa_path, k, paternal_kmers, maternal_kmers, path_delimiter);
    }
        // Min = 34 bits, max = 40 bits
    else if (k > 16 and k <= 20){
        phase<uint8_t,5>(gfa_path, k, paternal_kmers, maternal_kmers, path_delimiter);
    }
        // Min = 42 bits, max = 48 bits
    else if (k > 20 and k <= 24){
        phase<uint16_t,3>(gfa_path, k, paternal_kmers, maternal_kmers, path_delimiter);
    }
        // Min = 50 bits, max = 56 bits
    else if (k > 24 and k <= 28){
        phase<uint8_t,7>(gfa_path, k, paternal_kmers, maternal_kmers, path_delimiter);
    }
        // Min = 58 bits, max = 64 bits
    else if (k > 28 and k <= 32){
        phase<uint64_t,1>(gfa_path, k, paternal_kmers, maternal_kmers, path_delimiter);
    }
        // Min = 66 bits, max = 80 bits
    else if (k > 32 and k <= 40){
        phase<uint16_t,5>(gfa_path, k, paternal_kmers, maternal_kmers, path_delimiter);
    }
        // Min = 82 bits, max = 96 bits
    else if (k > 40 and k <= 48){
        phase<uint32_t,3>(gfa_path, k, paternal_kmers, maternal_kmers, path_delimiter);
    }
        // Min = 98 bits, max = 128 bits
    else if (k > 48 and k <= 64){
        phase<uint64_t,2>(gfa_path, k, paternal_kmers, maternal_kmers, path_delimiter);
    }
}

}
