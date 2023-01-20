#ifndef GFASE_BRIDGES_HPP
#define GFASE_BRIDGES_HPP

#include "bdsg/hash_graph.hpp"

#include <vector>
#include <functional>

using bdsg::HandleGraph;
using handlegraph::handle_t;

using std::vector;
using std::function;
using std::pair;

namespace gfase {

// return the bridge nodes of the handle graph using Schmidt's (2013)
// algorithm. assumes that graph is connected. the returned bridges are
// all in the forward orientation.
vector<handle_t> bridge_nodes(const HandleGraph& graph);

// consolidate runs of non-branching bridges into multi-node unipath bridges.
// the returned unipaths will be ordered / oriented so that they represent a
// valid walk of the bridge. the bridges argument can be a subset of all bridges.
// assumes bridges are all provided in the forward orientation.
vector<vector<handle_t>> consolidate_bridges(const HandleGraph& graph,
                                             const vector<handle_t>& bridges);


// execute a function on each bridge component, delimited by a (sub)set of consolidated
// unipath bridges. the function is also given the set of bridges that delimit
// the bridge component, which are identified by their index in the input vector
// and the direction from the bridge into the component (true = left). if a bridge
// consists of multiple nodes, then only the innermost one is included in the subgraph.
void for_each_bridge_component(const HandleGraph& graph,
                               const vector<vector<handle_t>>& bridges,
                               const function<void(const HandleGraph&,
                                                   const vector<pair<size_t, bool>>&)>& f);


}


#endif //GFASE_BRIDGES_HPP
