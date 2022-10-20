#include "MultiContactGraph.hpp"
#include "IncrementalIdMap.hpp"

using gfase::NonBipartiteEdgeException;
using gfase::IncrementalIdMap;
using gfase::MultiContactGraph;
using gfase::alt_component_t;

#include <iostream>

using std::exception;
using std::cerr;


int main(){
    {
        MultiContactGraph g;

        IncrementalIdMap<string> id_map(true);
        id_map.insert("a");
        id_map.insert("b");
        id_map.insert("c");

        g.insert_node(0);
        g.insert_node(1);
        g.insert_node(2);

        g.try_insert_edge(1,0,10);
        g.try_insert_edge(2,0,20);
        g.try_insert_edge(1,2,12);

        path output_path("test_contact_io.txt");
        g.write_contact_map(output_path, id_map);

        MultiContactGraph g2(output_path, id_map);

        path output_path2("test_contact_io2.txt");
        g2.write_contact_map(output_path2, id_map);
    }

    return 0;
}
