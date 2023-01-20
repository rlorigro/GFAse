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

void run_test(string& data_file, set<string>& correct_names) {
    
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();
    path relative_gfa_path = data_file;
    path absolute_gfa_path = project_directory / relative_gfa_path;
    
    GfaReader reader(absolute_gfa_path);
    
    HashGraph graph;
    IncrementalIdMap<string> id_map;
    
    gfa_to_handle_graph(graph, id_map, absolute_gfa_path);
    
    set<string> bridge_names;
    for_each_connected_component_subgraph(graph, [&](const HandleGraph& component) {
        for (auto bridge : bridge_nodes(component)) {
            bridge_names.insert(id_map.get_name(graph.get_id(bridge)));
        }
    });
    
    if (bridge_names != correct_names) {
        throw runtime_error(("ERROR: incorrect bridges detected in GFA " + data_file).c_str());
    }
}

int main(){

    {
        string file = "data/simple_chain.gfa";
        set<string> correct_names{"i", "j", "k", "l", "m"};
        run_test(file, correct_names);
    }
    {
        string file = "data/test_gfa1.gfa";
        set<string> correct_names{"11", "13"};
        run_test(file, correct_names);
    }
    {
        string file = "data/connected_components.gfa";
        set<string> correct_names{"a", "d", "e", "f", "g", "h"};
        run_test(file, correct_names);
    }
    {
        string file = "data/chain_test.gfa";
        // by component
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
        run_test(file, correct_names);
    }

    cerr << "All tests successful!" << endl;
    
    return 0;
}
