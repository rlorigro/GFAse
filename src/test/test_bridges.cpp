#include "Bridges.hpp"
#include "IncrementalIdMap.hpp"
#include "gfa_to_handle.hpp"
#include "handle_to_gfa.hpp"
#include "graph_utility.hpp"
#include "Filesystem.hpp"

#include "bdsg/hash_graph.hpp"

#include <string>
#include <set>
#include <iostream>
#include <algorithm>

using namespace gfase;

using ghc::filesystem::path;
using bdsg::HashGraph;
using handlegraph::handle_t;

using std::sort;
using std::string;
using std::cerr;
using std::endl;

void run_test(string& data_file,
              set<string>& correct_bridge_names,
              set<set<string>>& correct_bridge_components) {
    
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();
    path relative_gfa_path = data_file;
    path absolute_gfa_path = project_directory / relative_gfa_path;
    
    GfaReader reader(absolute_gfa_path);
    
    HashGraph graph;
    IncrementalIdMap<string> id_map;
    Overlaps overlaps;
    
    gfa_to_handle_graph(graph, id_map, overlaps, absolute_gfa_path);
    
    vector<handle_t> bridges;
    set<string> bridge_names;
    for_each_connected_component_subgraph(graph, [&](const HandleGraph& component) {
        for (auto bridge : bridge_nodes(component)) {
            bridges.push_back(bridge);
            bridge_names.insert(id_map.get_name(graph.get_id(bridge)));
        }
    });
    
    if (bridge_names != correct_bridge_names) {
        throw runtime_error(("ERROR: incorrect bridges detected in GFA " + data_file).c_str());
    }
    
    auto consolidated_bridges = consolidate_bridges(graph, bridges);
    
    set<set<string>> bridge_components;
    for_each_bridge_component(graph, consolidated_bridges,
                              [&](const HandleGraph& graph,
                                  const vector<pair<size_t, bool>>& incident_bridges) {
        set<string> component;
        
        graph.for_each_handle([&](const handle_t& handle) {
            component.insert(id_map.get_name(graph.get_id(handle)));
        });
        
        bridge_components.insert(component);
    });
    
    if (bridge_components != correct_bridge_components) {
        throw runtime_error(("ERROR: incorrect bridge components detected in GFA " + data_file).c_str());
    }
}

int main(){

    {
        string file = "data/simple_chain.gfa";
        set<string> correct_names{"i", "j", "k", "l", "m"};
        set<set<string>> correct_components{
            {"i"},
            {"m"},
            {"g", "h", "l", "m"},
            {"e", "f", "k", "l"},
            {"c", "d", "j", "k"},
            {"a", "b", "i", "j"}
        };
        run_test(file, correct_names, correct_components);
    }
    {
        string file = "data/test_gfa1.gfa";
        set<string> correct_names{"11", "13"};
        set<set<string>> correct_components{
            {"11", "12", "13"},
            {"11"},
            {"13"}
        };
        run_test(file, correct_names, correct_components);
    }
    {
        string file = "data/connected_components.gfa";
        set<string> correct_names{"a", "d", "e", "f", "g", "h"};
        set<set<string>> correct_components{
            {"a", "b", "c", "d"},
            {"a"},
            {"d"},
            {"g"},
            {"h"},
            {"e"},
            {"f"}
        };
        run_test(file, correct_names, correct_components);
    }
    {
        string file = "data/chain_test.gfa";
        set<string> correct_names{
            "a", "d", "g_a", "g_b", "j", "m", "n", "o", "p", "s", "t", "tip_a", "tip_d", "tip_e", "tip_f",
            "a2", "d2", "j2", "m2_0", "m2_1", "q2_0", "q2_1", "x2",
            "diploid_a0", "diploid_a1", "triploid_e0", "triploid_e1", "triploid_e2",
            "unique",
            "unphased_mat",
            "unphased_pat",
            "tree_a", "tree_e",
            "a3", "d3", "g3", "j3",
            "a4", "g4", "j4"
        };
        set<set<string>> correct_components{
            {"s"}, {"t"}, {"p", "q", "r", "s", "t"}, {"n"}, {"n", "o", "p"}, {"j", "k", "l", "m"}, {"g_b", "h", "i", "j"}, {"g_a", "g_b", "tip_a"}, {"tip_a", "tip_b", "tip_c", "tip_d"}, {"tip_d", "tip_e", "tip_f"}, {"tip_f"}, {"tip_e"}, {"d", "e", "f", "g_a"}, {"a", "b", "c", "d"}, {"a"},
            {"triploid_e1"}, {"triploid_e0"}, {"triploid_e2"}, {"diploid_a0", "diploid_a1", "diploid_b0", "diploid_b1", "diploid_c0", "diploid_c1", "diploid_d0", "diploid_d1", "diploid_e0", "diploid_e1", "triploid_e0", "triploid_e1", "triploid_e2"}, {"diploid_a1"}, {"diploid_a0"},
            {"q2_0"}, {"m2_0", "o2_0", "p2_0", "q2_0"}, {"k2_0", "k2_1", "l2_0", "l2_1", "m2_0", "m2_1", "x2"}, {"m2_1", "o2_1", "p2_1", "q2_1"}, {"q2_1"}, {"d2", "e2", "f2", "g2", "h2", "i2", "j2"}, {"a2", "b2", "c2", "d2"}, {"a2"},
            {"tree_e"}, {"tree_a", "tree_b0", "tree_b1", "tree_c0", "tree_c1", "tree_c2", "tree_c3", "tree_d0", "tree_d1", "tree_e"}, {"tree_a"},
            {"a3"}, {"a3", "b3", "c3", "d3"}, {"d3", "e3", "f3", "g3"}, {"g3", "h3", "i3", "j3"}, {"j3", "k3", "l3", "m3"},
            {"a4"}, {"a4", "b4", "c4", "e4", "f4", "g4"}, {"g4", "h4", "i4", "j4"}, {"j4", "k4", "l4", "m4"},
            {"unique"},
            {"unphased_mat"},
            {"unphased_pat"}
        };
        run_test(file, correct_names, correct_components);
    }

    cerr << "All tests successful!" << endl;
    
    return 0;
}
