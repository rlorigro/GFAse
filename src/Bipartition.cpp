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

        size_t subgraph_index = subgraphs.size();

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
            edge_t meta_edge = metagraph.edge_handle(h_a, h_b);
            meta_edge_to_edges[meta_edge].emplace_back(e);

            if (partition_a == partition_b){
                throw runtime_error("ERROR: direct edge between two subgraphs with the same partition membership");
            }
        }
    });
}


void Bipartition::for_each_subgraph(const function<void(const HandleGraph& subgraph, size_t subgraph_index, bool partition)>& f) const{
    for (size_t i=0; i<subgraphs.size(); i++){
        f(subgraphs[i], i, subgraph_partitions[i]);
    }
}


size_t Bipartition::get_subgraph_size(size_t subgraph_index) const{
    return subgraphs.at(subgraph_index).get_node_count();
}


size_t Bipartition::get_subgraph_index(nid_t meta_node) const{
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


void Bipartition::for_each_boundary_node_in_subgraph(size_t subgraph_index, bool left, const function<void(const handle_t& h)>& f){
    auto n = nid_t(subgraph_index);
    auto h = metagraph.get_handle(n);

    metagraph.follow_edges(h, left, [&](const handle_t& h_other){
        edge_t e;
        if (left) {
            e = metagraph.edge_handle(h_other, h);
        }
        else{
            e = metagraph.edge_handle(h, h_other);
        }

        for (auto& parent_edge: meta_edge_to_edges.at(e)){
            auto first_index = node_to_subgraph.at(graph.get_id(parent_edge.first));
            auto second_index = node_to_subgraph.at(graph.get_id(parent_edge.second));

            if (first_index == subgraph_index and second_index != subgraph_index){
                f(parent_edge.first);
            }
            else if (first_index != subgraph_index and second_index == subgraph_index){
                f(parent_edge.second);
            }
            else{
                throw runtime_error("ERROR: self-edge crosses subgraph boundary: " + to_string(first_index) + "->" + to_string(second_index));
            }
        }
    });

}


void Bipartition::for_each_handle_in_subgraph(size_t subgraph_index, const function<void(const handle_t& h)>& f){
    subgraphs.at(subgraph_index).for_each_handle(f);
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
}


void Bipartition::write_meta_graph_csv(ostream& file) const{
    file << "name" << ',' << "color" << ',' << "node_count" << ',' << '\n';

    for_each_subgraph([&](const HandleGraph& subgraph, size_t subgraph_index, bool partition){
        file << subgraph_index << ',' << (partition ? "#A2AFBE" : "#0D60BC") << ',' << subgraph.get_node_count() << ',' << '\n';
    });
}


}
