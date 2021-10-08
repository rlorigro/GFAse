#include "Bipartition.hpp"
#include "graph_utility.hpp"

namespace gfase {


Bipartition::Bipartition(PathHandleGraph& graph, IncrementalIdMap<string>& id_map):
    graph(graph),
    id_map(id_map)
{}


void Bipartition::partition(const unordered_set<nid_t>& node_subset){
    unordered_set<nid_t> all_nodes;

    graph.for_each_handle([&](const handle_t& h) {
        auto n = graph.get_id(h);
        all_nodes.emplace(n);
    });

    while (not all_nodes.empty()) {
        auto start_node = *all_nodes.begin();

        bool partition = 0;
        if (node_subset.find(start_node) == node_subset.end()){
            partition = 1;
        }

        size_t subgraph_index = subgraphs.size();

        // Allocate new elements in the vectors for this component
        subgraph_partitions.emplace_back(partition);
        subgraphs.emplace_back();

        // Enter all the nodes into a subgraph
        for_node_in_bfs(graph, start_node, node_subset, [&](const handle_t& h) {
            subgraphs.back().add_node(h);

            auto id = graph.get_id(h);
            node_to_subgraph.emplace(id, subgraph_index);
        });


        // TODO: figure out a way to save the edges and then create meta edges later
        // Duplicate all the edges
        for_edge_in_bfs(graph, start_node, node_subset, [&](const handle_t& handle_a, const handle_t& handle_b){},
            [&](const handle_t& handle_a, const handle_t& handle_b){
                auto id_a = graph.get_id(handle_a);
                auto id_b = graph.get_id(handle_b);

//                if (node_subset.find(id_a) == node_subset.end()){
//                    ;
//                }
//                else if (node_subset.find(id_b) == node_subset.end()){
//                    ;
//                }
//                else{
//                    throw runtime_error("ERROR: unexpected membership in do_not_visit for edge: " + name_a + " -> " + name_b);
//                }
            });
    }
}


}
