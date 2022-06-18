#include "IncrementalIdMap.hpp"
#include "ContactGraph.hpp"

using gfase::Node;
using gfase::ContactGraph;
using gfase::IncrementalIdMap;
using gfase::random_phase_search;

#include <iostream>

using std::cerr;

void test_mutability(){
    IncrementalIdMap<string> id_map(true);
    ContactGraph g;

    size_t n_nodes = 10;

    for (size_t i=0; i<n_nodes; i++) {
        string a = "n" + to_string(i);

        auto id = id_map.insert(a);
        g.insert_node(int32_t(id));
    }

    for (size_t i=0; i<n_nodes; i++){
        g.try_insert_edge(int32_t(i), int32_t((i+3) % n_nodes));
    }

    g.for_each_node([&](int32_t id, const Node& n){
        cerr << id << '\n';
        cerr << '\t' << "partition: " << int(n.partition) << '\n';
        cerr << '\t' << "neighbors: ";

        g.for_each_node_neighbor(id, [&](int32_t id_other, const Node& n_other){
            cerr << id_other << ' ';
        });

        cerr << '\n';
    });

    cerr << "Edges before editing:" << '\n';
    g.for_each_edge([&](const pair<int32_t,int32_t> e, int32_t weight){
        cerr << e.first << ',' << e.second << ' ' << int(weight) << '\n';
    });

    for (size_t i=5; i<n_nodes; i++){
        g.remove_edge(int32_t(i), int32_t((i+3) % n_nodes));
    }

    cerr << "Edges after editing:" << '\n';
    g.for_each_edge([&](const pair<int32_t,int32_t> e, int32_t weight){
        cerr << e.first << ',' << e.second << ' ' << int(weight) << '\n';
    });

    g.for_each_node([&](int32_t id, const Node& n){
        cerr << id << '\n';
        cerr << n << '\n';
        cerr << '\n';
    });

    for (size_t i=5; i<n_nodes; i++){
        g.try_insert_edge(int32_t(i), int32_t((i+3) % n_nodes));
    }

    cerr << "Edges after un-editing:" << '\n';
    g.for_each_edge([&](const pair<int32_t,int32_t> e, int32_t weight){
        cerr << e.first << ',' << e.second << ' ' << int(weight) << '\n';
    });

    g.for_each_node([&](int32_t id, const Node& n){
        cerr << id << '\n';
        cerr << n << '\n';
        cerr << '\n';
    });

    for (size_t i=5; i<n_nodes; i++){
        g.remove_node(int32_t(i));
    }

    cerr << "After removing nodes" << '\n';
    g.for_each_edge([&](const pair<int32_t,int32_t> e, int32_t weight){
        cerr << e.first << ',' << e.second << ' ' << int(weight) << '\n';
    });

    g.for_each_node([&](int32_t id, const Node& n){
        cerr << id << '\n';
        cerr << n << '\n';
        cerr << '\n';
    });

    g.for_each_edge([&](const pair<int32_t,int32_t> e, int32_t weight){
        g.increment_edge_weight(e.first, e.second, 777);
    });

    cerr << "After incrementing weights:" << '\n';
    g.for_each_edge([&](const pair<int32_t,int32_t> e, int32_t weight){
        cerr << e.first << ',' << e.second << ' ' << int(weight) << '\n';
    });
}


void test_optimization(){
    ContactGraph g;

    ///       (1)     (4)    (7)    (10)
    ///      /   \   /   \  /   \  /   \
    ///   (0)     (3)    (6)    (9)     (12)
    ///      \   /  \   /   \   /  \   /
    ///       (2)    (5)    (8)    (11)

    vector<int8_t> partitions = {0,0,-1,0,1,0,0,-1,-1,0,1,0,1};

    for (size_t i=0; i<13; i++){
        g.insert_node(int32_t(i),partitions[i]);
    }

    // Non-phaseable edges
    g.try_insert_edge(0,1,3);
    g.try_insert_edge(0,2,3);
    g.try_insert_edge(3,1,3);
    g.try_insert_edge(3,2,3);

    g.try_insert_edge(3,4,3);
    g.try_insert_edge(3,5,3);
    g.try_insert_edge(6,4,3);
    g.try_insert_edge(6,5,3);

    g.try_insert_edge(6,7,3);
    g.try_insert_edge(6,8,3);
    g.try_insert_edge(9,7,3);
    g.try_insert_edge(9,8,3);

    g.try_insert_edge(9,10,3);
    g.try_insert_edge(9,11,3);
    g.try_insert_edge(12,10,3);
    g.try_insert_edge(12,11,3);

    // Consistent edges
    g.try_insert_edge(1,4,6);
    g.try_insert_edge(1,7,3);
    g.try_insert_edge(1,10,1);

    g.try_insert_edge(4,7,6);
    g.try_insert_edge(4,10,3);

    g.try_insert_edge(7,10,6);

    // Consistent edges
    g.try_insert_edge(2,5,6);
    g.try_insert_edge(2,8,3);
    g.try_insert_edge(2,11,1);

    g.try_insert_edge(5,8,6);
    g.try_insert_edge(5,11,3);

    g.try_insert_edge(8,11,6);

    // Inconsistent edges
    g.try_insert_edge(1,5,1);
    g.try_insert_edge(4,8,1);
    g.try_insert_edge(7,11,1);

    g.try_insert_edge(2,4,1);
    g.try_insert_edge(5,7,1);
    g.try_insert_edge(8,10,1);

    // Intra-bubble edges
    g.try_insert_edge(1,2,1);
    g.try_insert_edge(4,5,1);
    g.try_insert_edge(7,8,1);
    g.try_insert_edge(10,11,1);

    cerr << "Before optimization:" << '\n';
    g.for_each_node([&](int32_t id, const Node& n){
        cerr << id << '\n';
        cerr << n << '\n';
        cerr << '\n';
    });

    vector <pair <int32_t,int8_t> > best_partitions;
    atomic<int64_t> best_score = std::numeric_limits<int64_t>::min();
    atomic<size_t> job_index = 0;
    mutex phase_mutex;
    size_t m_iterations = 100;

    random_phase_search(g, best_partitions, best_score, job_index, phase_mutex, m_iterations);

    cerr << "After optimization:" << '\n';
    g.for_each_node([&](int32_t id, const Node& n){
        cerr << id << '\n';
        cerr << n << '\n';
        cerr << '\n';
    });

    for (auto& [n,p]: best_partitions){
        cerr << n << ',' << int(p) << '\n';
    }

    // TODO: set/get partitions still broken!! graph internal state doesn't match "best_partitions"

}


int main(){
//    test_mutability();
    test_optimization();

    return 0;
}


