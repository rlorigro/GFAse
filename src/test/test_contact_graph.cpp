#include "IncrementalIdMap.hpp"
#include "ContactGraph.hpp"

using gfase::Node;
using gfase::ContactGraph;
using gfase::IncrementalIdMap;

#include <iostream>

using std::cerr;


int main(){
    IncrementalIdMap<string> id_map(true);
    ContactGraph g;

    size_t n_nodes = 10;

    for (size_t i=0; i<n_nodes; i++) {
        string a = "n" + to_string(i);

        auto id = id_map.insert(a);
        g.insert_node(int32_t(id));
    }

    for (size_t i=0; i<n_nodes; i++){
        g.insert_edge(int32_t(i), int32_t((i+3) % n_nodes));
    }

    g.for_each_node([&](int32_t id, const Node& n){
        cerr << id << '\n';
        cerr << '\t' << int(n.partition) << '\n';
        cerr << '\t';

        g.for_each_node_neighbor(id, [&](int32_t id_other, const Node& n_other){
            cerr << id_other << ',';
        });

        cerr << '\n';
    });

    g.for_each_edge([&](const pair<int32_t,int32_t> e, int32_t weight){
        cerr << e.first << ',' << e.second << ' ' << int(weight) << '\n';
    });

    return 0;
}


