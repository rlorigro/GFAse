#ifndef GFASE_GRAPH_UTILITY_HPP
#define GFASE_GRAPH_UTILITY_HPP

#include "IncrementalIdMap.hpp"
#include "SubgraphOverlay.hpp"
#include "gfa_to_handle.hpp"
#include "handle_to_gfa.hpp"
#include "Filesystem.hpp"
#include "GfaReader.hpp"
#include "Overlaps.hpp"
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

pair<string, size_t> parse_path_string(string path_name, char delimiter);

pair<nid_t,nid_t> translate_id(
        const HandleGraph& source_graph,
        const IncrementalIdMap<string>& source_id_map,
        HandleGraph& destination_graph,
        IncrementalIdMap<string>& destination_id_map,
        handle_t source_handle);

tuple<nid_t,nid_t,bool> try_translate_id(
        const HandleGraph& source_graph,
        const IncrementalIdMap<string>& source_id_map,
        HandleGraph& destination_graph,
        IncrementalIdMap<string>& destination_id_map,
        handle_t source_handle);

void for_node_in_bfs(
        const HandleGraph& graph,
        nid_t start_node,
        const function<bool(const nid_t& id)>& pass_criteria,
        const function<void(const handle_t& h)>& f);

void for_node_in_bfs(
        const HandleGraph& graph,
        nid_t start_node,
        const unordered_set<nid_t>& do_not_visit,
        const function<void(const handle_t&)>& f);

void for_node_in_bfs(const HandleGraph& graph, nid_t start_node, const function<void(const handle_t&)>& f);

void for_edge_in_bfs(
        const HandleGraph& graph,
        nid_t start_node,
        const unordered_set<nid_t>& do_not_visit,
        const function<void(const handle_t& handle_a, const handle_t& handle_b)>& f_pass,
        const function<void(const handle_t& handle_a, const handle_t& handle_b)>& f_fail);

void for_edge_in_bfs(
        const HandleGraph& graph,
        nid_t start_node,
        const function<void(const handle_t& handle_a, const handle_t& handle_b)>& f);

void for_each_connected_component(HandleGraph& graph, const function<void(unordered_set<nid_t>& connected_component)>& f);

void for_each_connected_component_subgraph(HandleGraph& graph, const function<void(const HandleGraph& subgraph)>& f);

void split_connected_components(
        MutablePathDeletableHandleGraph& graph,
        IncrementalIdMap<string>& id_map,
        Overlaps& overlaps,
        vector<HashGraph>& graphs,
        vector<IncrementalIdMap<string> >& id_maps,
        vector<Overlaps>& comp_overlaps,
        bool delete_visited_components = false);

void split_connected_components(
        MutablePathDeletableHandleGraph& graph,
        IncrementalIdMap<string>& id_map,
        Overlaps& overlaps,
        vector<HashGraph>& graphs,
        vector<IncrementalIdMap<string> >& id_maps,
        vector<Overlaps>& comp_overlaps,
        vector <vector <pair <string, string> > >& in_edges,
        vector <vector <pair <string, string> > >& out_edges,
        const unordered_set<nid_t>& do_not_visit,
        bool delete_visited_components);

void split_connected_components(
        MutablePathDeletableHandleGraph& graph,
        IncrementalIdMap<string>& id_map,
        vector<HashGraph>& graphs,
        bool delete_visited_components);

void write_connected_components_to_gfas(
        const MutablePathDeletableHandleGraph& graph,
        const IncrementalIdMap<string>& id_map,
        path output_directory);

void run_command(string& argument_string);

void plot_graph(const HandleGraph& graph, string filename_prefix);

void print_graph_paths(const PathHandleGraph& graph, const IncrementalIdMap<string>& id_map);

// Find any nodes that are adjacent to the beginning and end of a path, as long as they are the only adjacent node
pair<handle_t, bool> find_singleton_adjacent_handle(const PathHandleGraph& graph, const handle_t& h, bool left);

void find_diploid_paths(const PathHandleGraph& graph, unordered_set<string>& diploid_path_names);

void find_diploid_paths(
        const PathHandleGraph& graph,
        const set<string>& subset,
        unordered_set<string>& diploid_path_names,
        char path_delimiter);

void find_diploid_paths(
        const PathHandleGraph& graph,
        unordered_map<string, string>& diploid_path_names,
        unordered_set<string>& haploid_path_names);

void extend_paths(
        MutablePathMutableHandleGraph& graph,
        vector<pair<path_handle_t, handle_t> >& to_be_prepended,
        vector<pair<path_handle_t, handle_t> >& to_be_appended);

void un_extend_paths(
        MutablePathMutableHandleGraph& graph,
        const vector <pair<path_handle_t, handle_t> >& to_be_prepended,
        const vector <pair<path_handle_t, handle_t> >& to_be_appended);

void unzip(MutablePathDeletableHandleGraph& graph, IncrementalIdMap<string>& id_map, const Overlaps& overlaps, bool keep_paths=false, bool delete_islands=true);

void for_each_tip(const HandleGraph& graph, const function<void(const handle_t& h, bool is_left, bool is_right)>& f);

void write_paths_to_csv(const PathHandleGraph& graph, const IncrementalIdMap<string>& id_map, ofstream& file);

}

#endif //GFASE_GRAPH_UTILITY_HPP
