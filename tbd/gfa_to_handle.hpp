#ifndef GFASE_GFA_TO_HANDLE_HPP
#define GFASE_GFA_TO_HANDLE_HPP

/**
 * \file gfa_to_handle.hpp
 *
 * Defines algorithms for copying data from GFA files into handle graphs
 */

#include <iostream>
#include <cctype>
#include <string>

#include "gfakluge.hpp"
#include "bdsg/packed_graph.hpp"
#include "handlegraph/handle_graph.hpp"
#include "IncrementalIdMap.hpp"

using namespace std;
using std::runtime_error;
using handlegraph::nid_t;
using bdsg::MutableHandleGraph;
using bdsg::MutablePathMutableHandleGraph;
using gfase::IncrementalIdMap;

namespace gfase {

/// This exception will be thrown if the GFA data is not acceptable.
class GFAFormatError : public runtime_error {
public:
    using runtime_error::runtime_error;
};

/// Read a GFA file for a blunt-ended graph into a HandleGraph. Give "-" as a filename for stdin.
///
/// Optionally tries read the GFA from disk without creating an in-memory representation (defaults to
/// in-memory algorithm if reading from stdin).
///
/// Also optionally provides a hint about the node ID range to the handle graph implementation before
/// constructing it (defaults to no hint if reading from stdin).
///
/// Throws GFAFormatError if the GFA file is not acceptable, and
/// std::ios_base::failure if an IO operation fails. Throws invalid_argument if
/// otherwise misused.
void gfa_to_handle_graph(const string& filename,
                         MutableHandleGraph& graph,
                         IncrementalIdMap<string>& id_map,
                         bool try_from_disk = true,
                         bool try_id_increment_hint = false);

/// Same as gfa_to_handle_graph but also adds path elements from the GFA to the graph
void gfa_to_path_handle_graph(const string& filename,
                              MutablePathMutableHandleGraph& graph,
                              IncrementalIdMap<string>& id_map,
                              bool try_from_disk = true,
                              bool try_id_increment_hint = false);

/// Same as above but operating on a stream. Assumed to be non-seekable; all conversion happens in memory.
/// Always streaming. Doesn't support ID increment hints.
void gfa_to_path_handle_graph_in_memory(istream& in,
                                        MutablePathMutableHandleGraph& graph,
                                        IncrementalIdMap<string>& id_map);


}

#endif //GFASE_GFA_TO_HANDLE_HPP
