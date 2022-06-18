#include "ContactGraph.hpp"

#include <ostream>

using std::ostream;
using std::cerr;
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


void ContactGraph::insert_edge(int32_t a, int32_t b, int32_t weight){
    edge_weights.emplace(edge(a,b),weight);
    nodes.at(a).neighbors.emplace(b);
    nodes.at(b).neighbors.emplace(a);
}


void ContactGraph::try_insert_edge(int32_t a, int32_t b){
    if (nodes.count(a) == 0 or nodes.count(b) == 0){
        throw runtime_error("ERROR: cannot add edge with nonexistent node id: (" + to_string(a) + "," + to_string(b) +")");
    }

    auto e = edge(a,b);
    auto result = edge_weights.find(e);

    if (result == edge_weights.end()) {
        insert_edge(a,b,0);
    }
}


void ContactGraph::try_insert_edge(int32_t a, int32_t b, int32_t weight){
    if (nodes.count(a) == 0 or nodes.count(b) == 0){
        throw runtime_error("ERROR: cannot add edge with nonexistent node id: (" + to_string(a) + "," + to_string(b) +")");
    }

    auto e = edge(a,b);
    auto result = edge_weights.find(e);

    if (result == edge_weights.end()) {
        insert_edge(a,b,weight);
    }
}


void ContactGraph::increment_edge_weight(int32_t a, int32_t b, int32_t value){
    if (nodes.count(a) == 0 or nodes.count(b) == 0){
        throw runtime_error("ERROR: cannot add edge with nonexistent node id: (" + to_string(a) + "," + to_string(b) +")");
    }

    auto e = edge(a,b);
    auto result = edge_weights.find(e);

    if (result != edge_weights.end()) {
        result->second += value;
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


void ContactGraph::for_each_node_neighbor(int32_t id, const function<void(int32_t id_other, const Node& n)>& f) const{
    for (auto& id_other: nodes.at(id).neighbors){
        f(id_other, nodes.at(id_other));
    }
}


void ContactGraph::for_each_node(const function<void(int32_t id, const Node& n)>& f) const{
    for (int32_t id=0; id<nodes.size(); id++){
        f(id, nodes.at(id));
    }
}


void ContactGraph::for_each_edge(const function<void(const pair<int32_t,int32_t> edge, int32_t weight)>& f) const{
    for (auto& [e, weight]: edge_weights){
        f(e, weight);
    }
}


void ContactGraph::insert_node(int32_t id, int8_t partition){
    if (partition < -1 or partition > 1){
        throw runtime_error("ERROR: can't assign partition index outside of {-1,0,1}");
    }

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


void ContactGraph::try_insert_node(int32_t id, int8_t partition){
    if (nodes.count(id) == 0) {
        nodes.emplace(id, partition);
    }
}


void ContactGraph::set_partition(int32_t id, int8_t partition) {
    nodes.at(id).partition = partition;
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


size_t ContactGraph::size(){
    return nodes.size();
}


int64_t ContactGraph::compute_consistency_score(int32_t id) const{
    int64_t score = 0;

    auto n = nodes.at(id);

    for_each_node_neighbor(id, [&](int32_t id_other, const Node& n_other){
        score += n.partition * n_other.partition * edge_weights.at(edge(id, id_other));
    });

    return score;
}


int64_t ContactGraph::compute_total_consistency_score() const{
    int64_t score = 0;

    for_each_edge([&](const pair<int32_t,int32_t> edge, int32_t weight){
        score += nodes.at(edge.first).partition * nodes.at(edge.second).partition * weight;
    });

    return score;
}


void ContactGraph::get_partitions(vector <pair <int32_t,int8_t> >& partitions) const{
    partitions.clear();
    partitions.resize(nodes.size());

    size_t i=0;
    for (auto& [n, node]: nodes){
        partitions.at(i) = {n,node.partition};
        i++;
    }
}


void ContactGraph::set_partitions(const vector <pair <int32_t,int8_t> >& partitions){
    for (auto& [n, p]: partitions){
        set_partition(n, p);
    }
}


ostream& operator<<(ostream& o, const Node& n){
    o << '\t' << "partition: " << int(n.partition) << '\n';
    o << '\t' << "neighbors: ";

    for (auto& id: n.neighbors){
        o << id << ' ';
    }

    return o;
}


void random_phase_search(
        ContactGraph contact_graph,
        vector <pair <int32_t,int8_t> >& best_partitions,
        atomic<int64_t>& best_score,
        atomic<size_t>& job_index,
        mutex& phase_mutex,
        size_t m_iterations
){

    size_t m = job_index.fetch_add(1);

    // True random number
    std::random_device rd;

    // Pseudorandom generator with true random seed
    std::mt19937 rng(rd());
    std::uniform_int_distribution<int> uniform_distribution(0,int(contact_graph.size()-1));

    int64_t total_score;

    vector<size_t> order(contact_graph.size());
    for (size_t i=0; i<order.size(); i++){
        order[i] = i;
    }

    while (m < m_iterations) {
        // Randomly perturb
        for (size_t i=0; i < ((contact_graph.size()/5) + 1); i++) {
            auto partition = int8_t((uniform_distribution(rng) % 3) - 1);
            contact_graph.set_partition(uniform_distribution(rng), partition);
        }

        for (size_t i=0; i < contact_graph.size(); i++) {
            auto n = uniform_distribution(rng);

            int64_t max_score = std::numeric_limits<int64_t>::min();
            int8_t p_max = 0;

            for (int8_t p=-1; p<=1; p++) {
                contact_graph.set_partition(n, p);
                auto score = contact_graph.compute_consistency_score(n);

                if (score > max_score) {
                    max_score = score;
                    p_max = p;
                }
            }

            contact_graph.set_partition(n, p_max);
        }

        total_score = contact_graph.compute_total_consistency_score();

        phase_mutex.lock();
        if (total_score > best_score) {
            best_score = total_score;
            contact_graph.get_partitions(best_partitions);
        }
        else {
            contact_graph.set_partitions(best_partitions);
        }

        cerr << m << ' ' << best_score << ' ' << total_score << ' ' << '\n' << std::flush ;
        phase_mutex.unlock();

        m = job_index.fetch_add(1);
    }

    contact_graph.set_partitions(best_partitions);
}



}