#include "ContactGraph.hpp"

using std::min;
using std::max;


namespace gfase{


Node::Node(int8_t partition):
        neighbors(),
        partition(partition)
{}


ContactGraph::ContactGraph(const contact_map_t& contact_map, const IncrementalIdMap<string>& id_map) {
    string s;

    for (const auto& [a,sub_map]: contact_map){
        try_insert_node(a);

        for (const auto& [b,count]: sub_map){
            try_insert_node(b);
            insert_edge(a,b,count);
        }
    }
}


pair<int32_t,int32_t> edge(int32_t a, int32_t b){
    return {min(a,b), max(a,b)};
}


void ContactGraph::insert_edge(int32_t a, int32_t b){
    if (nodes.count(a) == 0 or nodes.count(b) == 0){
        throw runtime_error("ERROR: cannot add edge with nonexistent node id: (" + to_string(a) + "," + to_string(b) +")");
    }

    auto e = edge(a,b);
    auto result = edge_weights.find(e);

    if (result == edge_weights.end()) {
        edge_weights.emplace(e,0);
        nodes.at(a).neighbors.emplace(b);
        nodes.at(b).neighbors.emplace(a);
    }
}


void ContactGraph::insert_edge(int32_t a, int32_t b, int32_t weight){
    if (nodes.count(a) == 0 or nodes.count(b) == 0){
        throw runtime_error("ERROR: cannot add edge with nonexistent node id: (" + to_string(a) + "," + to_string(b) +")");
    }

    auto e = edge(a,b);
    auto result = edge_weights.find(e);

    if (result == edge_weights.end()) {
        edge_weights.emplace(e,weight);
        nodes.at(a).neighbors.emplace(b);
        nodes.at(b).neighbors.emplace(a);
    }
}


void ContactGraph::remove_edge(int32_t a, int32_t b){
    auto e = edge(a,b);
    auto result = edge_weights.find(e);

    if (result != edge_weights.end()) {
        edge_weights.erase(result);
        nodes.at(a).neighbors.erase(b);
        nodes.at(b).neighbors.erase(a);
    }
}


void ContactGraph::for_each_node_neighbor(int32_t node_id, const function<void(int32_t id, const Node& n)>& f){
    for (auto& id: nodes.at(node_id).neighbors){
        f(id, nodes.at(node_id));
    }
}


void ContactGraph::for_each_node(const function<void(int32_t id, const Node& n)>& f){
    for (int32_t id=0; id<nodes.size(); id++){
        f(id, nodes.at(id));
    }
}


void ContactGraph::for_each_edge(const function<void(const pair<int32_t,int32_t> edge, int32_t weight)>& f){
    for (auto& [e, weight]: edge_weights){
        f(e, weight);
    }
}


void ContactGraph::insert_node(int32_t id, int8_t partition){
    nodes.emplace(id, partition);
}


void ContactGraph::insert_node(int32_t id){
    nodes.emplace(id, 0);
}


void ContactGraph::try_insert_node(int32_t id){
    if (nodes.count(id) == 0) {
        nodes.emplace(id, 0);
    }
}


void ContactGraph::remove_node(int32_t id){
    vector <pair <int32_t, int32_t> > to_be_removed;

    for_each_node_neighbor(id, [&](int32_t id_other, const Node& n){
        to_be_removed.emplace_back(edge(id, id_other));
    });

    for (auto& e: to_be_removed){
        remove_edge(e.first, e.second);
    }

    nodes.erase(id);
}


}