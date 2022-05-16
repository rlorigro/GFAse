#include "Bipartition.hpp"
#include "graph_utility.hpp"

namespace gfase {


Bipartition::Bipartition(PathHandleGraph& graph, IncrementalIdMap<string>& id_map, unordered_set<nid_t>& node_subset):
    graph(graph),
    id_map(id_map),
    node_subset(node_subset)
{}


bool Bipartition::get_partition_of_node(const nid_t& id) const{
    return node_subset.find(id) == node_subset.end();
}


void Bipartition::partition(){
    unordered_set<nid_t> all_nodes;
    unordered_set<edge_t> all_edges;

    graph.for_each_handle([&](const handle_t& h) {
        auto n = graph.get_id(h);
        all_nodes.emplace(n);
    });

    while (not all_nodes.empty()) {
        auto start_node = *all_nodes.begin();
        bool partition = get_partition_of_node(start_node);

        // Skip 0 and start at 1 because the world is already on fire anyway
        max_id++;

        size_t subgraph_index = max_id;

        // Allocate new elements in the hashmap for this component
        subgraph_partitions.emplace(subgraph_index, partition);
        subgraphs.emplace(subgraph_index, &graph);

        metagraph.create_handle("", nid_t(subgraph_index));

        // Enter all the nodes into a subgraph
        for_node_in_bfs(graph, start_node,
                [&](const nid_t& id){
                    return get_partition_of_node(id) == partition;
                },
                [&](const handle_t& h) {
                    subgraphs[subgraph_index].add_node(h);

                    auto id = graph.get_id(h);
                    node_to_subgraph.emplace(id, subgraph_index);

                    all_nodes.erase(id);
                });
    }

    graph.for_each_edge([&](const edge_t& e){
        auto e_normalized = graph.edge_handle(e.first, e.second);

        auto id_a = graph.get_id(e_normalized.first);
        auto id_b = graph.get_id(e_normalized.second);

        auto subgraph_index_a = node_to_subgraph.at(id_a);
        auto subgraph_index_b = node_to_subgraph.at(id_b);

        auto partition_a = get_partition_of_node(id_a);
        auto partition_b = get_partition_of_node(id_b);

        if (subgraph_index_a != subgraph_index_b){
            auto h_a = metagraph.get_handle(nid_t(subgraph_index_a), graph.get_is_reverse(e_normalized.first));
            auto h_b = metagraph.get_handle(nid_t(subgraph_index_b), graph.get_is_reverse(e_normalized.second));

            auto meta_edge = metagraph.edge_handle(h_a, h_b);
            metagraph.create_edge(meta_edge.first, meta_edge.second);

            // The metagraph edge can be used to find all individual edges that linked subgraphs in the parent graph
            meta_edge_to_edges[meta_edge].emplace(e_normalized);

            if (partition_a == partition_b){
                throw runtime_error("ERROR: direct edge between two subgraphs with the same partition membership");
            }
        }
    });
}


void Bipartition::merge_subgraphs(size_t subgraph_index_a, size_t subgraph_index_b){

    if (subgraph_partitions.at(subgraph_index_b) != subgraph_partitions.at(subgraph_index_a)){
        throw runtime_error("ERROR: cannot merge subgraphs of differing partitions");
    }

    subgraphs.at(subgraph_index_b).for_each_handle([&](const handle_t& h){
        // Add all the b nodes to a
        subgraphs.at(subgraph_index_a).add_node(h);

        // Overwrite mappings from parent to subgraph
        node_to_subgraph[graph.get_id(h)] = subgraph_index_a;
    });

    auto metagraph_handle_b = metagraph.get_handle(nid_t(subgraph_index_b));
    auto metagraph_handle_a = metagraph.get_handle(nid_t(subgraph_index_a));

    unordered_set<edge_t> to_be_deleted;

    // Duplicate the metagraph boundary edges to a and remove the existing ones from b (LEFT)
    metagraph.follow_edges(metagraph_handle_b, true, [&](const handle_t & h){
        edge_t e_a = metagraph.edge_handle(h, metagraph_handle_a);
        edge_t e_b = metagraph.edge_handle(h, metagraph_handle_b);

        for (auto& parent_edge: meta_edge_to_edges.at(e_b)){
            meta_edge_to_edges[e_a].emplace(graph.edge_handle(parent_edge.first,parent_edge.second));
        }

        to_be_deleted.emplace(e_b);
    });

    // Duplicate the metagraph boundary edges to a and remove the existing ones from b (RIGHT)
    metagraph.follow_edges(metagraph_handle_b, false, [&](const handle_t & h){
        edge_t e_a = metagraph.edge_handle(metagraph_handle_a, h);
        edge_t e_b = metagraph.edge_handle(metagraph_handle_b, h);

        for (auto& parent_edge: meta_edge_to_edges.at(e_b)){
            meta_edge_to_edges[e_a].emplace(graph.edge_handle(parent_edge.first,parent_edge.second));
        }

        to_be_deleted.emplace(e_b);
    });

    for (auto& e: to_be_deleted){
        meta_edge_to_edges.erase(e);
    }

    subgraphs.erase(subgraph_index_b);
    metagraph.destroy_handle(metagraph_handle_b);
    subgraph_partitions.erase(subgraph_index_b);
}


void Bipartition::for_each_subgraph(const function<void(const HandleGraph& subgraph, size_t subgraph_index, bool partition)>& f) const{
    for (auto& [sugraph_index, subgraph]: subgraphs){
        f(subgraph, sugraph_index, subgraph_partitions.at(sugraph_index));
    }
}


size_t Bipartition::get_subgraph_size(size_t subgraph_index) const{
    return subgraphs.at(subgraph_index).get_node_count();
}


size_t Bipartition::get_subgraph_index_of_parent_node(nid_t meta_node) const{
    return node_to_subgraph.at(meta_node);
}


nid_t Bipartition::get_id(handle_t meta_handle) const{
    return metagraph.get_id(meta_handle);
}


void Bipartition::follow_subgraph_edges(size_t subgraph_index, bool go_left, const function<void(const handle_t& h)>& f){
    auto n = nid_t(subgraph_index);
    auto h = metagraph.get_handle(n);

    metagraph.follow_edges(h, go_left, [&](const handle_t& h_other){
        edge_t e;
        if (go_left) {
            e = metagraph.edge_handle(h_other, h);
        }
        else{
            e = metagraph.edge_handle(h, h_other);
        }

        for (auto& parent_edge: meta_edge_to_edges.at(e)){
            auto first_index = node_to_subgraph.at(graph.get_id(parent_edge.first));
            auto second_index = node_to_subgraph.at(graph.get_id(parent_edge.second));

            if (first_index == subgraph_index and second_index != subgraph_index){
                f(parent_edge.second);
            }
            else if (first_index != subgraph_index and second_index == subgraph_index){
                f(parent_edge.first);
            }
            else{
                throw runtime_error("ERROR: self-edge crosses subgraph boundary: " + to_string(first_index) + "->" + to_string(second_index));
            }
        }
    });

}


void Bipartition::for_each_boundary_node_in_subgraph(size_t subgraph_index, bool left, const function<void(const handle_t& h)>& f) const{
    auto n = nid_t(subgraph_index);
    auto h = metagraph.get_handle(n, false);

    unordered_set <handle_t> boundary_nodes;

//    cerr << "Looking for boundary nodes of subgraph: " << subgraph_index << '\n';

    for (auto& l: {true,false}){
//        cerr << (l ? "checking left" : "checking right") << '\n';

        metagraph.follow_edges(h, l, [&](const handle_t& h_other){
            edge_t e;
            if (l) {
                e = metagraph.edge_handle(h_other, h);
            }
            else{
                e = metagraph.edge_handle(h, h_other);
            }

//            cerr << "\tmeta-edge: "
//                 << graph.get_id(e.first) << (graph.get_is_reverse(e.first) ? '-' : '+') << " -> "
//                 << graph.get_id(e.second) << (graph.get_is_reverse(e.second) ? '-' : '+') << '\n';

            for (auto parent_edge: meta_edge_to_edges.at(e)){
//                cerr << "\t\tparent edge: "
//                     << id_map.get_name(graph.get_id(parent_edge.first)) << (graph.get_is_reverse(parent_edge.first) ? '-' : '+') << " -> "
//                     << id_map.get_name(graph.get_id(parent_edge.second)) << (graph.get_is_reverse(parent_edge.second) ? '-' : '+') << '\n';

                auto first_index = node_to_subgraph.at(graph.get_id(parent_edge.first));
                auto second_index = node_to_subgraph.at(graph.get_id(parent_edge.second));

                if (first_index == subgraph_index and second_index != subgraph_index){
                    // Verify that the edge has the handle-of-interest in the F orientation
                    if (graph.get_is_reverse(parent_edge.first)){
                        if (left){
//                            cerr << "\t\t\tfound: " << id_map.get_name(graph.get_id(parent_edge.first)) << '\n';
                            boundary_nodes.emplace(graph.flip(parent_edge.first));
                        }
                    }
                    else if (not left){
                        boundary_nodes.emplace(parent_edge.first);
                    }
                }
                else if (first_index != subgraph_index and second_index == subgraph_index){
                    // Verify that the edge has the handle-of-interest in the F orientation
                    if (graph.get_is_reverse(parent_edge.second)){
                        if (not left){
//                            cerr << "\t\t\tfound: " << id_map.get_name(graph.get_id(parent_edge.second)) << '\n';
                            boundary_nodes.emplace(graph.flip(parent_edge.second));
                        }
                    }
                    else if (left){
                        boundary_nodes.emplace(parent_edge.second);
                    }
                }
                else{
                    throw runtime_error("ERROR: self-edge crosses subgraph boundary: " + to_string(first_index) + "->" + to_string(second_index));
                }
            }
        });
    }

    for (auto& item: boundary_nodes){
        f(item);
    }
}


void Bipartition::for_each_handle_in_subgraph(size_t subgraph_index, const function<void(const handle_t& h)>& f){
    subgraphs.at(subgraph_index).for_each_handle(f);
}


nid_t Bipartition::get_id_of_parent_handle(const handle_t& h){
    return graph.get_id(h);
}


nid_t Bipartition::get_id_of_parent_handle(const string& name){
    return id_map.get_id(name);
}


string Bipartition::get_name_of_parent_node(nid_t id){
    return id_map.get_name(id);
}


size_t Bipartition::size() const{
    return subgraphs.size();
}


bool Bipartition::get_partition_of_subgraph(const size_t subgraph_index) const{
    return subgraph_partitions.at(subgraph_index);
}


void Bipartition::print() const{
    for_each_subgraph([&](const HandleGraph& subgraph, size_t subgraph_index, bool partition){
        cerr << "Subgraph: " << subgraph_index << '\n';
        cerr << "\tPartition: " << int(partition) << '\n';
        cerr << "\tNodes: " << subgraph.get_node_count() << '\n';

        subgraph.for_each_handle([&](const handle_t& h){
            // Debug printing
            auto id = subgraph.get_id(h);
            string name = id_map.get_name(id);

            cerr << '\t' << '\t' << name << '\n';
        });
    });
}


void Bipartition::write_parent_graph_csv(ostream& file) const{
    file << "name" << ',' << "color" << ',' << "subgraph_index" << ',' << '\n';

    for_each_subgraph([&](const HandleGraph& subgraph, size_t subgraph_index, bool partition){
        subgraph.for_each_handle([&](const handle_t& h){
            auto id = subgraph.get_id(h);
            string name = id_map.get_name(id);

            file << name << ',' << (partition ? "#A2AFBE" : "#0D60BC") << ',' << subgraph_index << ',' << '\n';
        });
    });

    file << std::flush;
}


void Bipartition::write_meta_graph_csv(ostream& file) const{
    file << "name" << ',' << "color" << ',' << "node_count" << ',' << '\n';

    for_each_subgraph([&](const HandleGraph& subgraph, size_t subgraph_index, bool partition){
        file << subgraph_index << ',' << (partition ? "#A2AFBE" : "#0D60BC") << ',' << subgraph.get_node_count() << ',' << '\n';
    });

    file << std::flush;
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

                // Make sure it isn't possible to visit more than 2 phased bubble sides from this unphased node
                unordered_set<size_t> left_first_degree_neighbors;
                unordered_set<size_t> right_first_degree_neighbors;

                ploidy_bipartition.follow_subgraph_edges(subgraph_index, true, [&](const handle_t& h){
                    auto id = ploidy_bipartition.get_id(h);
                    auto adjacent_subgraph_index = ploidy_bipartition.get_subgraph_index_of_parent_node(id);

                    if (ploidy_bipartition.get_partition_of_subgraph(adjacent_subgraph_index) == 0) {
                        left_first_degree_neighbors.emplace(adjacent_subgraph_index);
                    }

                    ploidy_bipartition.follow_subgraph_edges(adjacent_subgraph_index, false, [&](const handle_t& h2){
                        auto id2 = ploidy_bipartition.get_id(h2);
                        auto adjacent_subgraph_index2 = ploidy_bipartition.get_subgraph_index_of_parent_node(id2);
                        second_degree_neighbors.emplace(adjacent_subgraph_index2);
                    });
                });

                ploidy_bipartition.follow_subgraph_edges(subgraph_index, false, [&](const handle_t& h){
                    auto id = ploidy_bipartition.get_id(h);
                    auto adjacent_subgraph_index = ploidy_bipartition.get_subgraph_index_of_parent_node(id);

                    if (ploidy_bipartition.get_partition_of_subgraph(adjacent_subgraph_index) == 0) {
                        right_first_degree_neighbors.emplace(adjacent_subgraph_index);
                    }

                    ploidy_bipartition.follow_subgraph_edges(adjacent_subgraph_index, true, [&](const handle_t& h2){
                        auto id2 = ploidy_bipartition.get_id(h2);
                        auto adjacent_subgraph_index2 = ploidy_bipartition.get_subgraph_index_of_parent_node(id2);
                        second_degree_neighbors.emplace(adjacent_subgraph_index2);
                    });
                });

                // If there are no second degree neighbors, this unphased subgraph passes
                if (second_degree_neighbors.size() == 1 and right_first_degree_neighbors.size() < 3 and left_first_degree_neighbors.size() < 3){
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


void for_element_in_bubble_chain(
        const Bipartition& chain_bipartition,
        const HandleGraph& graph,
        const IncrementalIdMap<string>& id_map,
        const BubbleGraph& bubble_graph,
        const function<void(const vector<string>& node_names, size_t subgraph_index)>& f
){

    chain_bipartition.for_each_subgraph([&](const HandleGraph& subgraph, size_t subgraph_index, bool partition){
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
//            cerr << "queuing start node: " << id_map.get_name(graph.get_id(h)) << (graph.get_is_reverse(h) ? '-' : '+') << '\n';
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

        // Initialize things for this chain
        if (not next_nodes.empty() and subgraph.get_node_count() > 1){
            q.emplace(next_nodes);
        }

        // Iterate each bubble or bridge and update the queue with the next nodes
        while(not q.empty()){
            auto& nodes = q.front();

            if (nodes.size() == 1){
                auto node = *nodes.begin();
                auto id = subgraph.get_id(node);

                // Verify not diploid by name
                auto name = id_map.get_name(id);

                // Single node elements in a chain should not be flagged as diploid
                if (bubble_graph.node_is_bubble(int32_t(id))){
                    throw runtime_error("ERROR: non-bubble in chain is flagged as diploid: " + name);
                }
                else {
                    vector<string> n = {name};
                    f(n, subgraph_index);
                }
            }
            else if (nodes.size() == 2){
                // Verify nodes are diploid counterparts to one another
                auto& node_a = *nodes.begin();
                auto& node_b = *(++nodes.begin());

                auto id_a = subgraph.get_id(node_a);
                auto id_b = subgraph.get_id(node_b);

                auto name_a = id_map.get_name(id_a);
                auto name_b = id_map.get_name(id_b);

                if (not bubble_graph.node_is_bubble(int32_t(id_a))){
                    throw runtime_error("ERROR: non diploid node in bubble of bubble chain: " + name_a);
                }

                if (bubble_graph.get_other_side(int32_t(id_a)) != id_b){
                    throw runtime_error("ERROR: nodes in bubble are not labeled as diploid counterparts: " + name_a + "," + name_b);
                }

                vector<string> n = {name_a, name_b};
                f(n, subgraph_index);
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


}
