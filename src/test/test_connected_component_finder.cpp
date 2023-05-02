#include "IncrementalIdMap.hpp"
#include "gfa_to_handle.hpp"
#include "handle_to_gfa.hpp"
#include "graph_utility.hpp"
#include "Filesystem.hpp"

#include "bdsg/hash_graph.hpp"

#include <string>

using gfase::IncrementalIdMap;
using gfase::handle_graph_to_gfa;
using gfase::for_each_connected_component;
using gfase::split_connected_components;
using gfase::print_graph_paths;
using gfase::plot_graph;

using ghc::filesystem::path;
using bdsg::HashGraph;
using bdsg::MutablePathMutableHandleGraph;
using handlegraph::path_handle_t;
using handlegraph::step_handle_t;
using handlegraph::handle_t;

using std::string;
using std::cout;
using std::cerr;


int main(){
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();

    // Get test VCF path
    path relative_gfa_path = "data/connected_components.gfa";
    path absolute_gfa_path = project_directory / relative_gfa_path;

    HashGraph graph;
    IncrementalIdMap<string> id_map;
    Overlaps overlaps;

    gfa_to_handle_graph(graph, id_map, overlaps, absolute_gfa_path);

    unordered_set<string> cc1 = {"a", "b", "c"};
    unordered_set<string> cc2 = {"d", "e"};
    bool found_cc1 = false;
    bool found_cc2 = false;

    for_each_connected_component(graph, [&](unordered_set<nid_t>& connected_component){
        unordered_set<string> cc;

        for (auto& n: connected_component){
            auto name = id_map.get_name(n);
            cc.emplace(name);
        }

        if (cc == cc1){
            found_cc1 = true;
        }
        else if (cc == cc2){
            found_cc2 = true;
        }
        else{
            cerr << "BAD COMPONENT:" << '\n';
            for (auto& item: cc){
                cerr << item << '\n';
            }
            throw runtime_error("FAIL: connected component does not match any in truth set");
        }
    });

    if (not found_cc1){
        throw runtime_error("FAIL: cc1 not found");
    }

    if (not found_cc2){
        throw runtime_error("FAIL: cc2 not found");
    }

    vector<HashGraph> connected_component_graphs;
    vector <IncrementalIdMap<string> > connected_component_ids;
    vector<Overlaps> connected_component_overlaps;

    split_connected_components(graph, id_map, overlaps, connected_component_graphs, connected_component_ids, connected_component_overlaps);

    for (size_t i=0; i<connected_component_graphs.size(); i++) {
        plot_graph(connected_component_graphs[i], "component_" + to_string(i));
        print_graph_paths(connected_component_graphs[i], connected_component_ids[i]);
    }

    return 0;
}

