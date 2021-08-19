#include "IncrementalIdMap.hpp"
#include "gfa_to_handle.hpp"
#include "handle_to_gfa.hpp"
#include "GraphUtility.hpp"
#include "Filesystem.hpp"

#include "bdsg/hash_graph.hpp"

#include <string>

using gfase::IncrementalIdMap;
using gfase::handle_graph_to_gfa;
using gfase::for_node_in_bfs;
using gfase::for_edge_in_bfs;

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
    path relative_gfa_path = "data/simple_chain.gfa";
    path absolute_gfa_path = project_directory / relative_gfa_path;

    GfaReader reader(absolute_gfa_path);

    HashGraph graph;
    IncrementalIdMap<string> id_map;

    gfa_to_handle_graph(graph, id_map, absolute_gfa_path);

    set<string> bfs_node_names;
    for_node_in_bfs(graph, 1, [&](const handle_t& h){
        auto n = graph.get_id(h);
        auto name = id_map.get_name(n);

        bfs_node_names.emplace(name);
        cerr << "Iterating node: " << name << '\n';
    });

    set<string> all_node_names;
    reader.for_each_sequence([&](string& name, string& sequence){
        all_node_names.emplace(name);
    });

    cerr << "All nodes in GFA\n";
    for (auto& item: all_node_names){
        cerr << item << '\n';
    }

    cerr << "All nodes in BFS\n";
    for (auto& item: bfs_node_names){
        cerr << item << '\n';
    }

    if (not (all_node_names == bfs_node_names)){
        vector<string> bfs_only;
        set_difference(bfs_node_names.begin(), bfs_node_names.end(), all_node_names.begin(), all_node_names.end(), std::inserter(bfs_only, bfs_only.begin()));

        vector<string> gfa_only;
        set_difference(all_node_names.begin(), all_node_names.end(), bfs_node_names.begin(), bfs_node_names.end(), std::inserter(gfa_only, gfa_only.begin()));

        cerr << "Nodes found only in BFS\n";
        for (auto& item: bfs_only){
            cerr << item << '\n';
        }

        cerr << "Nodes found only in GFA\n";
        for (auto& item: gfa_only){
            cerr << item << '\n';
        }

        throw runtime_error("FAIL: bfs node names do not match GfaReader node names");
    }

    graph.for_each_handle([&](const handle_t& h){
        auto id = graph.get_id(h);
        auto name = id_map.get_name(id);

        cerr << id << ' ' << name << '\n';
    });

//    set <pair <string,string> > bfs_edge_names;
    for_edge_in_bfs(graph, 1, [&](const handle_t& h1, const handle_t& h2){
        auto n1 = graph.get_id(h1);
        auto name1 = id_map.get_name(n1);

        auto n2 = graph.get_id(h2);
        auto name2 = id_map.get_name(n2);

//        bfs_node_names.emplace(name1, name2);
        cerr << "Iterating edge: " << name1 << (graph.get_is_reverse(h1) ? '-' : '+') << ' ' << name2 << (graph.get_is_reverse(h2) ? '-' : '+') << '\n';
    });


    return 0;
}
