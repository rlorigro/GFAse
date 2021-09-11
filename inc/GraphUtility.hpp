#ifndef GFASE_GRAPHUTILITY_HPP
#define GFASE_GRAPHUTILITY_HPP

#include "IncrementalIdMap.hpp"
#include "gfa_to_handle.hpp"
#include "handle_to_gfa.hpp"
#include "Filesystem.hpp"
#include "GfaReader.hpp"
#include "CLI11.hpp"

#include "bdsg/hash_graph.hpp"

#include <functional>
#include <string>
#include <queue>

using gfase::IncrementalIdMap;
using gfase::handle_graph_to_gfa;

using ghc::filesystem::path;
using bdsg::MutablePathMutableHandleGraph;
using bdsg::MutablePathDeletableHandleGraph;
using bdsg::PathHandleGraph;
using bdsg::path_handle_t;
using bdsg::step_handle_t;
using bdsg::handle_t;
using bdsg::HashGraph;

using std::string;
using std::queue;
using std::cout;
using std::cerr;

namespace gfase{

void for_node_in_bfs(const HandleGraph& graph, nid_t start_node, const function<void(const handle_t&)>& f);

void for_edge_in_bfs(const HandleGraph& graph, nid_t start_node, const function<void(const handle_t& handle_a, const handle_t& handle_b)>& f);

void for_each_connected_component(HandleGraph& graph, const function<void(unordered_set<nid_t>& connected_component)>& f);

void split_connected_components(
        MutablePathDeletableHandleGraph& graph,
        IncrementalIdMap<string>& id_map,
        vector<HashGraph>& graphs,
        vector<IncrementalIdMap<string> >& id_maps,
        bool delete_visited_components = false);

void write_connected_components_to_gfas(
        const MutablePathDeletableHandleGraph& graph,
        const IncrementalIdMap<string>& id_map,
        path output_directory);

void run_command(string& argument_string);

void plot_graph(const HandleGraph& graph, string filename_prefix);

void print_graph_paths(const PathHandleGraph& graph, const IncrementalIdMap<string>& id_map);


}

#endif //GFASE_GRAPHUTILITY_HPP
