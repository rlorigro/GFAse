#ifndef GFASE_HANDLE_TO_GFA_HPP
#define GFASE_HANDLE_TO_GFA_HPP

#include "handlegraph/path_handle_graph.hpp"
#include "handlegraph/handle_graph.hpp"
#include "IncrementalIdMap.hpp"
#include <fstream>

using bdsg::PathHandleGraph;
using bdsg::HandleGraph;
using bdsg::path_handle_t;
using bdsg::handle_t;
using bdsg::edge_t;
using std::runtime_error;
using std::ostream;
using std::string;


namespace gfase {


char get_reversal_character(const HandleGraph& graph, const handle_t& node);

void write_node_to_gfa(const HandleGraph& graph, const handle_t& node, ostream& output_file);

void write_node_to_gfa(const HandleGraph& graph, const IncrementalIdMap<string>& id_map, const handle_t& node, ostream& output_file);

void write_edge_to_gfa(const HandleGraph& graph, const edge_t& edge, ostream& output_file);

void write_edge_to_gfa(const HandleGraph& graph, const IncrementalIdMap<string>& id_map, const edge_t& edge, ostream& output_file);

void write_path_to_gfa(
        const PathHandleGraph& graph,
        const IncrementalIdMap<string>& id_map,
        const path_handle_t& path,
        ostream& output_file);

void handle_graph_to_gfa(const HandleGraph& graph, ostream& output_gfa);

void handle_graph_to_gfa(const PathHandleGraph& graph, const IncrementalIdMap<string>& id_map, ostream& output_gfa);

}

#endif //GFASE_HANDLE_TO_GFA_HPP
