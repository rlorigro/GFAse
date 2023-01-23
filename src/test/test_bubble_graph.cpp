#include "IncrementalIdMap.hpp"
#include "BubbleGraph.hpp"

using gfase::Bubble;
using gfase::BubbleGraph;
using gfase::IncrementalIdMap;

#include <iostream>

using std::cerr;


int main(){
    IncrementalIdMap<string> id_map;
    BubbleGraph bg;

    for (size_t i=0; i<100; i++) {
        auto a = 2*i;
        auto b = 2*i + 1;

        string name_a = to_string(a);
        string name_b = to_string(b);

        auto id_a = id_map.insert(name_a);
        auto id_b = id_map.insert(name_b);

        auto result = bg.try_add_bubble(int32_t(id_a), int32_t(id_b));

        if (i % 3 == 0){
            bg.flip(result);
        }
    }

    path output_path = "test_phases.csv";
    bg.write_bandage_csv(output_path, id_map);

    IncrementalIdMap<string> id_map_2;
    BubbleGraph bg2(output_path, id_map_2);

    vector<bool> p;
    vector<bool> p2;

    // When the CSV gets loaded, the state of the bubbles is lost in terms of whether it is "flipped", but
    // the phases assigned to each ID should be identical
    for (size_t i=0; i<100; i++) {
        auto a = 2*i;
        auto b = 2*i + 1;

        string name_a = to_string(a);
        string name_b = to_string(b);

        auto old_id_a = id_map.get_id(name_a);
        auto new_id_a = id_map_2.get_id(name_a);

        auto old_bubble_id = bg.find_bubble_id_of_node(int32_t(old_id_a));
        auto new_bubble_id = bg2.find_bubble_id_of_node(int32_t(new_id_a));

        Bubble<int32_t> old_bubble;
        Bubble<int32_t> new_bubble;
        bg.at(new_bubble_id, new_bubble);
        bg.at(old_bubble_id, old_bubble);

        cerr << new_bubble.get(0) << ',' << new_bubble.get(1) << " " << old_bubble.get(0) << ',' << old_bubble.get(1) << '\n';

        if (new_bubble.get(0) != old_bubble.get(0)){
            throw runtime_error("ERROR: old bubble does not match new bubble");
        }
    }

    return 0;
}
