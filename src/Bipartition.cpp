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

    vector<edge_t> to_be_deleted;

    // Duplicate the metagraph boundary edges to a and remove the existing ones from b (LEFT)
    metagraph.follow_edges(metagraph_handle_b, true, [&](const handle_t & h){
        edge_t e_a = {h, metagraph_handle_a};
        edge_t e_b = {h, metagraph_handle_b};

        for (auto& parent_edge: meta_edge_to_edges.at(e_b)){
            meta_edge_to_edges[e_a].emplace_back(parent_edge);
        }

        to_be_deleted.emplace_back(e_b);
    });

    // Duplicate the metagraph boundary edges to a and remove the existing ones from b (RIGHT)
    metagraph.follow_edges(metagraph_handle_b, false, [&](const handle_t & h){
        edge_t e_a = {metagraph_handle_a, h};
        edge_t e_b = {metagraph_handle_b, h};

        for (auto& parent_edge: meta_edge_to_edges.at(e_b)){
            meta_edge_to_edges[e_a].emplace_back(parent_edge);
        }

        to_be_deleted.emplace_back(e_b);
    });

    for (auto& e: to_be_deleted){
        meta_edge_to_edges.erase(e);
    }

    subgraphs.erase(subgraph_index_b);
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


}
