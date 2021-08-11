#ifndef GFASE_GFA_TO_HANDLE_HPP
#define GFASE_GFA_TO_HANDLE_HPP

/**
 * \file gfa_to_handle.hpp
 *
 * Defines algorithms for copying data from GFA files into handle graphs
 */

#include "bdsg/hash_graph.hpp"
#include "handlegraph/handle_graph.hpp"
#include "IncrementalIdMap.hpp"
#include "GfaReader.hpp"

#include <cctype>
#include <string>
#include <iostream>

using namespace std;
using std::runtime_error;
using handlegraph::nid_t;
using handlegraph::step_handle_t;
using bdsg::MutableHandleGraph;
using bdsg::MutablePathMutableHandleGraph;
using gfase::IncrementalIdMap;

namespace gfase {

nid_t parse_gfa_sequence_id(const string& s, IncrementalIdMap<string>& id_map);

void gfa_to_handle_graph(
        MutablePathMutableHandleGraph& graph,
        IncrementalIdMap<string>& id_map,
        unordered_map<nid_t,step_handle_t>& node_to_path_step,
        path gfa_file_path);



}

#endif //GFASE_GFA_TO_HANDLE_HPP
