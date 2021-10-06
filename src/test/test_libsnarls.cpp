#include "IncrementalIdMap.hpp"
#include "gfa_to_handle.hpp"
#include "handle_to_gfa.hpp"
#include "graph_utility.hpp"
#include "Filesystem.hpp"

#include "snarls/handle_graph_snarl_finder.hpp"
#include "snarls/integrated_snarl_finder.hpp"
#include "snarls/snarl_finder.hpp"
#include "bdsg/hash_graph.hpp"

#include <string>

using gfase::IncrementalIdMap;
using gfase::handle_graph_to_gfa;
using gfase::for_node_in_bfs;
using gfase::for_edge_in_bfs;
using gfase::plot_graph;

using ghc::filesystem::path;
using bdsg::HashGraph;
using bdsg::MutablePathMutableHandleGraph;
using handlegraph::path_handle_t;
using handlegraph::step_handle_t;
using handlegraph::handle_t;
using snarls::IntegratedSnarlFinder;
using snarls::HandleGraphSnarlFinder;
using snarls::SnarlFinder;

using std::string;
using std::cout;
using std::cerr;


int main(){
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();

    // Get test VCF path
    path relative_gfa_path = "data/simple_chain.gfa";
    path absolute_gfa_path = project_directory / relative_gfa_path;

    HashGraph graph;
    IncrementalIdMap<string> id_map;

    gfa_to_handle_graph(graph, id_map, absolute_gfa_path);

    plot_graph(graph, "start_graph");

    IntegratedSnarlFinder snarl_finder(graph);
    vector <pair<nid_t, bool> > stack;

    // Stolen from Xian Chang
    // https://github.com/xchang1/vg/blob/snarl-traversal-bug/src/subcommand/stats_main.cpp#L1037-L1075
    snarl_finder.traverse_decomposition(
            [&](handle_t chain_start_handle) {
                nid_t id = graph.get_id(chain_start_handle);
                bool rev = graph.get_is_reverse(chain_start_handle);
                stack.emplace_back(id, rev);
                cerr << string(stack.size() - 1,'\t');
                cerr << "start chain " << id_map.get_name(id) << " " << rev << endl;
            },
            [&](handle_t chain_end_handle) {
                nid_t id = graph.get_id(chain_end_handle);
                bool rev = graph.get_is_reverse(chain_end_handle);
                pair<nid_t, bool> start = stack.back();
                cerr << string(stack.size() - 1,'\t');
                stack.pop_back();
                cerr << "end chain " << id_map.get_name(start.first) << " " << start.second << "->" << id_map.get_name(id) << " " << rev << endl;
            },
            [&](handle_t snarl_start_handle) {
                nid_t id = graph.get_id(snarl_start_handle);
                bool rev = graph.get_is_reverse(snarl_start_handle);
                stack.emplace_back(id, rev);
                cerr << string(stack.size() - 1,'\t');
                cerr << "start snarl " << id_map.get_name(id) << " " << rev << endl;
            },
            [&](handle_t snarl_end_handle) {
                nid_t id = graph.get_id(snarl_end_handle);
                bool rev = graph.get_is_reverse(snarl_end_handle);
                pair<nid_t, bool> start = stack.back();
                cerr << string(stack.size() - 1,'\t');
                stack.pop_back();
                cerr << "end snarl " << id_map.get_name(start.first) << " " << start.second << "->" << id_map.get_name(id) << " " << rev << endl;
            }
    );

    return 0;
}
