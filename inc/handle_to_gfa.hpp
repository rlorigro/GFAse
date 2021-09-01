#ifndef GFASE_HANDLE_TO_GFA_HPP
#define GFASE_HANDLE_TO_GFA_HPP

#include "handlegraph/path_handle_graph.hpp"
#include "handlegraph/handle_graph.hpp"
#include "IncrementalIdMap.hpp"
#include <fstream>

using handlegraph::PathHandleGraph;
using handlegraph::HandleGraph;
using handlegraph::handle_t;
using handlegraph::edge_t;
using std::runtime_error;
using std::ostream;
using std::string;


namespace gfase {


char get_reversal_character(const HandleGraph& graph, const handle_t& node);


void write_node_to_gfa(const HandleGraph& graph, const handle_t& node, ostream& output_file);


void write_edge_to_gfa(const HandleGraph& graph, const edge_t& edge, ostream& output_file);


/// With no consideration for directionality, just dump all the edges/nodes into GFA format
void handle_graph_to_gfa(const HandleGraph& graph, ostream& output_gfa);

void handle_graph_to_gfa(const PathHandleGraph& graph, const IncrementalIdMap<string>& id_map, ostream& output_gfa);

// TODO write this method to use the overlaps and id map to write the linkages/sequences in the canonical direction
// using the canonical names as well, wherever possible
void handle_graph_to_canonical_gfa(const HandleGraph& graph, const string& output_path);

}

#endif //GFASE_HANDLE_TO_GFA_HPP
