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
              const set<vector<pair<string, bool>>>& correct_phase_paths) {
    
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();
    path relative_gfa_path = data_file;
    path absolute_gfa_path = project_directory / relative_gfa_path;
    
    GfaReader reader(absolute_gfa_path);
    
    HashGraph graph;
    IncrementalIdMap<string> id_map;
    
    gfa_to_handle_graph(graph, id_map, absolute_gfa_path);
    
    cerr << "node ID translation:" << endl;
    graph.for_each_handle([&](const handle_t& handle) {
        cerr << graph.get_id(handle) << " -> " << id_map.get_name(graph.get_id(handle)) << endl;
    });
    
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
    for (auto node_set : {make_pair(phase_0_nodes, true), make_pair(phase_1_nodes, false)}) {
        for (auto node_name : node_set.first) {
            int32_t node_id = id_map.get_id(node_name);
            contact_graph.insert_node(node_id, node_set.second ? -1 : 1);
        }
    }
    for (auto alt_pair : alt_pairs) {
        contact_graph.add_alt(id_map.get_id(alt_pair.first), id_map.get_id(alt_pair.second));
    }
    
    set<vector<handle_t>> correct_phase_handle_paths;
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
        correct_phase_handle_paths.insert(phase_handle_path);
    }
    
    HamiltonianChainer chainer;
    chainer.generate_chain_paths(graph, contact_graph);
    
    set<vector<handle_t>> identified_phase_handle_paths;
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
        identified_phase_handle_paths.emplace(steps);
    });
    
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
        set<vector<pair<string, bool>>> correct_phase_paths{
            {{"i", false}, {"a", false}, {"j", false}, {"d", false}, {"k", false}, {"f", false}, {"l", false}},
            {{"i", false}, {"b", false}, {"j", false}, {"c", false}, {"k", false}, {"e", false}, {"l", false}}
        };
        run_test(file, phase_0_nodes, phase_1_nodes, alt_pairs, correct_phase_paths);
    }
    
    return 0;
}
