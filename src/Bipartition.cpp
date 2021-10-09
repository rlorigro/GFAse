#include "Bipartition.hpp"
#include "graph_utility.hpp"

namespace gfase {


Bipartition::Bipartition(PathHandleGraph& graph, IncrementalIdMap<string>& id_map, unordered_set<nid_t>& node_subset):
    graph(graph),
    id_map(id_map),
    node_subset(node_subset)
{}


bool Bipartition::get_partition_of_node(nid_t id){
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

        size_t subgraph_index = subgraphs.size();

        cerr << "New subgraph: " << subgraph_index << '\n';

        // Allocate new elements in the vectors for this component
        subgraph_partitions.emplace_back(partition);
        subgraphs.emplace_back(&graph);
        metagraph.create_handle("", nid_t(subgraph_index));

        // Enter all the nodes into a subgraph
        for_node_in_bfs(graph, start_node,
                [&](const nid_t& id){
                    return get_partition_of_node(id) == partition;
                },
                [&](const handle_t& h) {
                    subgraphs.back().add_node(h);

                    auto id = graph.get_id(h);
                    node_to_subgraph.emplace(id, subgraph_index);

                    cerr << '\t' << id_map.get_name(id) << '\n';

                    all_nodes.erase(id);
                });
    }

    graph.for_each_edge([&](const edge_t& e){
        auto id_a = graph.get_id(e.first);
        auto id_b = graph.get_id(e.second);

        auto subgraph_index_a = node_to_subgraph.at(id_a);
        auto subgraph_index_b = node_to_subgraph.at(id_b);

        auto partition_a = get_partition_of_node(id_a);
        auto partition_b = get_partition_of_node(id_b);

        if (subgraph_index_a != subgraph_index_b){
            auto h_a = metagraph.get_handle(nid_t(subgraph_index_a), graph.get_is_reverse(e.first));
            auto h_b = metagraph.get_handle(nid_t(subgraph_index_b), graph.get_is_reverse(e.second));

            metagraph.create_edge(h_a, h_b);

            // The metagraph edge can be used to find all individual edges that linked subgraphs in the parent graph
            edge_t meta_edge(h_a, h_b);
            meta_edge_to_edges[meta_edge].emplace_back(e);

            if (partition_a == partition_b){
                throw runtime_error("ERROR: direct edge between two subgraphs with the same partition membership");
            }
        }
    });
}


void Bipartition::for_each_subgraph(const function<void(const HandleGraph& subgraph, size_t subgraph_index, bool partition)>& f){
    for (size_t i=0; i<subgraphs.size(); i++){
        f(subgraphs[i], i, subgraph_partitions[i]);
    }
}


}
