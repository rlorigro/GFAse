#include "HamiltonianChainer.hpp"
#include "IncrementalIdMap.hpp"
#include "MultiContactGraph.hpp"
#include "gfa_to_handle.hpp"
#include "handle_to_gfa.hpp"
#include "graph_utility.hpp"
#include "Filesystem.hpp"

#include "bdsg/hash_graph.hpp"
#include "handlegraph/util.hpp"

#include <string>
#include <set>
#include <iostream>
#include <algorithm>
#include <vector>
#include <utility>
#include <cassert>

using namespace gfase;

using ghc::filesystem::path;
using bdsg::HashGraph;
using handlegraph::handle_t;
using handlegraph::nid_t;

using std::sort;
using std::string;
using std::cerr;
using std::endl;
using std::set;
using std::vector;
using std::pair;

void run_test(string& data_file,
              const set<string>& phase_0_nodes,
              const set<string>& phase_1_nodes,
              const set<pair<string, string>> alt_pairs,
              const vector<vector<pair<string, bool>>>& correct_phase_paths) {
    
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();
    path relative_gfa_path = data_file;
    path absolute_gfa_path = project_directory / relative_gfa_path;
    
    GfaReader reader(absolute_gfa_path);
    
    HashGraph graph;
    IncrementalIdMap<string> id_map;
    Overlaps overlaps;
    
    gfa_to_handle_graph(graph, id_map, overlaps, absolute_gfa_path);
    
//    cerr << "node ID translation:" << endl;
//    graph.for_each_handle([&](const handle_t& handle) {
//        cerr << graph.get_id(handle) << " -> " << id_map.get_name(graph.get_id(handle)) << endl;
//    });
    
    // get rid of pre-existing paths
    vector<path_handle_t> embedded_paths;
    graph.for_each_path_handle([&](const path_handle_t& path_handle) {
        embedded_paths.push_back(path_handle);
    });
    for (auto path_handle : embedded_paths) {
        graph.destroy_path(path_handle);
    }
    
    // fix the partition and alt relationships
    MultiContactGraph contact_graph;
    for (auto node_set : {phase_0_nodes, phase_1_nodes}) {
        for (auto node_name : node_set) {
            int32_t node_id = id_map.get_id(node_name);
            contact_graph.insert_node(node_id);
        }
    }
    for (auto alt_pair : alt_pairs) {
        contact_graph.add_alt(id_map.get_id(alt_pair.first), id_map.get_id(alt_pair.second));
    }
    vector<alt_component_t> alt_components;
    contact_graph.get_alt_components(alt_components);
    for (auto& alt_component : alt_components) {
        unordered_set<int> partitions_first, partitions_second;
        for (auto set_pair : {make_pair(alt_component.first, &partitions_first), make_pair(alt_component.second, &partitions_second)}) {
            for (auto node_id : set_pair.first) {
                int partition;
                if (phase_0_nodes.count(id_map.get_name(node_id))) {
                    partition = -1;
                }
                else if (phase_1_nodes.count(id_map.get_name(node_id))) {
                    partition = 1;
                }
                else {
                    partition = 0;
                }
                set_pair.second->insert(partition);
            }
        }
        assert(!partitions_first.empty() && !partitions_second.empty());
        assert(partitions_first.empty() || partitions_first.size() == 1);
        assert(partitions_second.empty() || partitions_second.size() == 1);
        if (!partitions_second.empty() && !partitions_first.empty()) {
            assert(*partitions_first.begin() == -(*partitions_second.begin()));
        }
        contact_graph.set_partition(*alt_component.first.begin(), *partitions_first.begin());
    }
    
    vector<vector<handle_t>> correct_phase_handle_paths;
    for (const auto& phase_path : correct_phase_paths) {
        vector<handle_t> phase_handle_path;
        if (id_map.get_id(phase_path.front().first) < id_map.get_id(phase_path.back().first)) {
            for (auto it = phase_path.begin(); it != phase_path.end(); ++it) {
                phase_handle_path.push_back(graph.get_handle(id_map.get_id(it->first), it->second));
            }
        }
        else {
            for (auto it = phase_path.rbegin(); it != phase_path.rend(); ++it) {
                phase_handle_path.push_back(graph.get_handle(id_map.get_id(it->first), !it->second));
            }
        }
        correct_phase_handle_paths.push_back(phase_handle_path);
    }
    sort(correct_phase_handle_paths.begin(), correct_phase_handle_paths.end());
    
    HamiltonianChainer chainer;
    chainer.generate_chain_paths(graph, id_map, contact_graph);
    
    vector<vector<handle_t>> identified_phase_handle_paths;
    graph.for_each_path_handle([&](const path_handle_t& path_handle) {
        vector<handle_t> steps;
        for (handle_t step : graph.scan_path(path_handle)) {
            steps.push_back(step);
        }
        if (graph.get_id(steps.back()) < graph.get_id(steps.front())) {
            vector<handle_t> rev_steps;
            for (auto it = steps.rbegin(); it != steps.rend(); ++it) {
                rev_steps.push_back(graph.flip(*it));
            }
            steps = rev_steps;
        }
        identified_phase_handle_paths.emplace_back(steps);
    });
    sort(identified_phase_handle_paths.begin(), identified_phase_handle_paths.end());
    
    if (identified_phase_handle_paths != correct_phase_handle_paths) {
        for (bool correct : {true, false}) {
            cerr << (correct ? "expected" : "obtained") << " phase paths:" << endl;
            for (auto phase_path : (correct ? correct_phase_handle_paths : identified_phase_handle_paths)) {
                for (auto handle : phase_path) {
                    cerr << " " << graph.get_id(handle) << (graph.get_is_reverse(handle) ? "-" : "+");
                }
                cerr << endl;
            }
        }
        throw runtime_error(("ERROR: incorrect phase paths detected in GFA " + data_file).c_str());
    }
}

int main(){

    {
        string file = "data/simple_chain_long_haploid.gfa";
        set<string> phase_0_nodes{"a", "d", "f"};
        set<string> phase_1_nodes{"b", "c", "e"};
        set<pair<string, string>> alt_pairs{{"a", "b"}, {"c", "d"}, {"e", "f"}};
        vector<vector<pair<string, bool>>> correct_phase_paths{
            {{"i", false}, {"a", false}, {"j", false}, {"d", false}, {"k", false}, {"f", false}, {"l", false}},
            {{"i", false}, {"b", false}, {"j", false}, {"c", false}, {"k", false}, {"e", false}, {"l", false}}
        };
        run_test(file, phase_0_nodes, phase_1_nodes, alt_pairs, correct_phase_paths);
    }
    {
        string file = "data/simple_chain_empty_node.gfa";
        set<string> phase_0_nodes{"a", "d", "f", "g"};
        set<string> phase_1_nodes{"b", "c", "e", "h"};
        set<pair<string, string>> alt_pairs{{"a", "b"}, {"c", "d"}, {"e", "f"}, {"g", "h"}};
        vector<vector<pair<string, bool>>> correct_phase_paths{
            {{"i", false}, {"a", false}, {"j", false}, {"d", false}, {"k", false}, {"f", false}, {"l", false}, {"g", false}, {"m", false}},
            {{"i", false}, {"b", false}, {"j", false}, {"c", false}, {"k", false}, {"e", false}, {"l", false}, {"h", false}, {"m", false}}
        };
        run_test(file, phase_0_nodes, phase_1_nodes, alt_pairs, correct_phase_paths);
    }
    {
        // induce a phase break by an unphased bubble
        string file = "data/simple_chain_empty_node.gfa";
        set<string> phase_0_nodes{"a", "d", "g"};
        set<string> phase_1_nodes{"b", "c", "h"};
        set<pair<string, string>> alt_pairs{{"a", "b"}, {"c", "d"}, {"g", "h"}};
        vector<vector<pair<string, bool>>> correct_phase_paths{
            {{"i", false}, {"a", false}, {"j", false}, {"d", false}, {"k", false}},
            {{"i", false}, {"b", false}, {"j", false}, {"c", false}, {"k", false}},
            {{"l", false}, {"g", false}, {"m", false}},
            {{"l", false}, {"h", false}, {"m", false}}
        };
        run_test(file, phase_0_nodes, phase_1_nodes, alt_pairs, correct_phase_paths);
    }
    {
        // include unphased, bridge-only components
        string file = "data/connected_components.gfa";
        set<string> phase_0_nodes{"b"};
        set<string> phase_1_nodes{"c"};
        set<pair<string, string>> alt_pairs{{"b", "c"}};
        vector<vector<pair<string, bool>>> correct_phase_paths{
            {{"a", false}, {"b", false}, {"d", false}},
            {{"a", false}, {"c", false}, {"d", false}}
        };
        run_test(file, phase_0_nodes, phase_1_nodes, alt_pairs, correct_phase_paths);
    }
    {
        string file = "data/big_test.gfa";
        set<string> phase_0_nodes{
            "b", "e",
            "h", "n", "q", "s",
            "tip_b", "tip_e",
            "b3", "f3", "h3", "k3",
            "a2", "b2", "f2", "h2",
            "k2_0", "o2_0",
            "k2_1", "o2_1"
        };
        set<string> phase_1_nodes{
            "c", "f",
            "i", "o", "r", "t",
            "tip_c", "tip_f",
            "c3", "e3", "i3", "l3",
            "c2", "e2", "i2",
            "l2_0", "p2_0",
            "l2_1", "p2_1"
        };
        set<pair<string, string>> alt_pairs{
            {"b", "c"}, {"e", "f"},
            {"h", "i"}, {"n", "o"}, {"q", "r"}, {"s", "t"},
            {"tip_b", "tip_c"}, {"tip_e", "tip_f"},
            {"c3", "b3"}, {"f3", "e3"}, {"h3", "i3"}, {"k3", "l3"},
            {"a2", "c2"}, {"b2", "c2"}, {"e2", "f2"}, {"h2", "i2"},
            {"k2_0", "l2_0"}, {"o2_0", "p2_0"},
            {"k2_1", "l2_1"}, {"o2_1", "p2_1"}
        };
        vector<vector<pair<string, bool>>> correct_phase_paths{
            {{"a", false}, {"b", false}, {"d", false}, {"e", false}, {"g_a", false}},
            {{"a", false}, {"c", false}, {"d", false}, {"f", false}, {"g_a", false}},
            {{"g_b", false}, {"h", false}, {"j", false}},
            {{"g_b", false}, {"i", false}, {"j", false}},
            {{"m", false}, {"n", false}, {"p", false}, {"q", false}, {"s", false}},
            {{"m", false}, {"o", false}, {"p", false}, {"r", false}, {"t", false}},
            {{"tip_a", false}, {"tip_b", false}, {"tip_d", false}, {"tip_e", false}},
            {{"tip_a", false}, {"tip_c", false}, {"tip_d", false}, {"tip_f", false}},
            {{"a3", false}, {"b3", false}, {"d3", false}, {"f3", false}, {"g3", false}, {"h3", false}, {"j3", false}, {"k3", false}}, // TODO: this misses m3 because it's not a bridge...
            {{"a3", false}, {"c3", false}, {"d3", false}, {"e3", false}, {"g3", false}, {"i3", false}, {"j3", false}, {"l3", false}},
            {{"a2", false}, {"b2", false}, {"d2", false}, {"f2", false}, {"g2", false}, {"h2", false}, {"j2", false}, {"x2", false}},
            {{"c2", false}, {"d2", false}, {"e2", false}, {"g2", false}, {"i2", false}, {"j2", false}, {"x2", false}},
            {{"l2_0", false}, {"m2_0", false}, {"p2_0", false}, {"q2_0", false}},
            {{"k2_0", false}, {"m2_0", false}, {"o2_0", false}, {"q2_0", false}},
            {{"l2_1", false}, {"m2_1", false}, {"p2_1", false}, {"q2_1", false}},
            {{"k2_1", false}, {"m2_1", false}, {"o2_1", false}, {"q2_1", false}}
        };
        run_test(file, phase_0_nodes, phase_1_nodes, alt_pairs, correct_phase_paths);
    }
    {
        // TODO: the current algorithm can't really do anything with the triploid component, but
        // maybe it could if we inserted a dummy node into non-star biclique adjacency components,
        // which could then act as a bridge
        string file = "data/chain_test.gfa";
        set<string> phase_0_nodes{
            "b", "e",
            "tip_c", "tip_f",
            "h", "k", "n", "q", "s",
            "b2", "f2", "h2",
            "k2_0", "o2_0",
            "k2_1", "o2_1",
            "b3", "f3", "h3", "k3",
            "tree_c2", "tree_c0",
            "c4", "f4", "i4", "k4"
        };
        set<string> phase_1_nodes{
            "c", "f",
            "tip_b", "tip_e",
            "i", "l", "o", "r", "t",
            "c2", "e2", "i2",
            "l2_0", "p2_0",
            "l2_1", "p2_1",
            "c3", "e3", "i3", "l3",
            "tree_c3", "tree_c1",
            "b4", "e4", "h4", "l4"
        };
        set<pair<string, string>> alt_pairs{
            {"b", "c"}, {"e", "f"},
            {"tip_b", "tip_c"}, {"tip_e", "tip_f"},
            {"h", "i"}, {"k", "l"}, {"n", "o"}, {"q", "r"}, {"s", "t"},
            {"b2", "c2"}, {"e2", "f2"}, {"h2", "i2"},
            {"k2_0", "l2_0"}, {"o2_0", "p2_0"},
            {"k2_1", "l2_1"}, {"o2_1", "p2_1"},
            {"c3", "b3"}, {"f3", "e3"}, {"h3", "i3"}, {"k3", "l3"},
            {"tree_c2", "tree_c3"}, {"tree_c0", "tree_c1"},
            {"c4", "b4"}, {"e4", "f4"}, {"h4", "i4"}, {"k4", "l4"}
        };
        vector<vector<pair<string, bool>>> correct_phase_paths{
            {{"a", false}, {"b", false}, {"d", false}, {"e", false}, {"g_a", false}},
            {{"a", false}, {"c", false}, {"d", false}, {"f", false}, {"g_a", false}},
            {{"tip_a", false}, {"tip_b", false}, {"tip_d", false}, {"tip_e", false}},
            {{"tip_a", false}, {"tip_c", false}, {"tip_d", false}, {"tip_f", false}},
            {{"g_b", false}, {"h", false}, {"j", false}, {"k", false}, {"m", false}},
            {{"n", false}, {"p", false}, {"q", false}, {"s", false}},
            {{"g_b", false}, {"i", false}, {"j", false}, {"l", false}, {"m", false}, {"o", false}, {"p", false}, {"r", false}, {"t", false}},
            {{"a2", false}, {"b2", false}, {"d2", false}, {"f2", false}},
            {{"h2", false}, {"j2", false}, {"x2", false}},
            {{"a2", false}, {"c2", false}, {"d2", false}, {"e2", false}},
            {{"i2", false}, {"j2", false}, {"x2", false}},
            {{"l2_0", false}, {"m2_0", false}, {"p2_0", false}, {"q2_0", false}},
            {{"k2_0", false}, {"m2_0", false}, {"o2_0", false}, {"q2_0", false}},
            {{"l2_1", false}, {"m2_1", false}, {"p2_1", false}, {"q2_1", false}},
            {{"k2_1", false}, {"m2_1", false}, {"o2_1", false}, {"q2_1", false}},
            {{"a3", false}, {"b3", false}, {"d3", false}, {"f3", false}, {"g3", false}, {"h3", false}, {"j3", false}, {"k3", false}},
            {{"a3", false}, {"c3", false}, {"d3", false}, {"e3", false}, {"g3", false}, {"i3", false}, {"j3", false}, {"l3", false}},
            {{"tree_b1", false}, {"tree_c3", false}, {"tree_d1", false}},
            {{"tree_b1", false}, {"tree_c2", false}, {"tree_d1", false}},
            {{"tree_b0", false}, {"tree_c0", false}, {"tree_d0", false}},
            {{"tree_b0", false}, {"tree_c1", false}, {"tree_d0", false}},
            {{"a4", false}, {"b4", false}, {"e4", true}, {"g4", false}, {"h4", false}, {"j4", true}, {"l4", false}},
            {{"a4", false}, {"c4", false}, {"f4", false}, {"g4", false}, {"i4", true}, {"j4", true}, {"k4", false}}
        };
        run_test(file, phase_0_nodes, phase_1_nodes, alt_pairs, correct_phase_paths);
    }
    
    cerr << "All tests successful!" << endl;
    
    return 0;
}
