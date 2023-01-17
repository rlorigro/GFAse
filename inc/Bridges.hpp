#ifndef GFASE_BRIDGES_HPP
#define GFASE_BRIDGES_HPP

#include "hash_graph.hpp"

#include <vector>

using bdsg::HandleGraph;
using handlegraph::handle_t;

using std::vector;

namespace gfase {

// return the bridge nodes of the handle graph using Schmidt's (2013)
// algorithm, assumes that graph is connected
vector<handle_t> bridge_nodes(const HandleGraph& graph);


}


#endif //GFASE_BRIDGES_HPP
