#include "IncrementalIdMap.hpp"
#include <iostream>

using gfase::IncrementalIdMap;
using std::cerr;

int main(){
    {
        IncrementalIdMap<string> id_map(true);

        id_map.insert("a");
        id_map.insert("b");
        id_map.insert("c");

        path output_path = "test_id_map.txt";
        id_map.write_to_csv(output_path);

        IncrementalIdMap<string> id_map_2(output_path);

        for (const auto&[name, id]: id_map_2.ids) {
            cerr << id << ',' << name << '\n';
        }
    }

    {
        IncrementalIdMap<string> id_map(false);

        id_map.insert("a");
        id_map.insert("b");
        id_map.insert("c");

        path output_path = "test_id_map.txt";
        id_map.write_to_csv(output_path);

        IncrementalIdMap<string> id_map_2(output_path);

        for (const auto&[name, id]: id_map_2.ids) {
            cerr << id << ',' << name << '\n';
        }
    }

    return 0;
}