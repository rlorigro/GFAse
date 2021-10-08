#ifndef GFASE_BIPARTITION_HPP
#define GFASE_BIPARTITION_HPP

#include "IncrementalIdMap.hpp"

#include "handlegraph/path_handle_graph.hpp"
#include "bdsg/overlays/packed_subgraph_overlay.hpp"
#include "bdsg/hash_graph.hpp"

#include <unordered_set>
#include <array>

using handlegraph::PathHandleGraph;
using handlegraph::edge_t;
using bdsg::PackedSubgraphOverlay;
using bdsg::HashGraph;

using std::unordered_set;
using std::array;


namespace gfase{

class Bipartition {
public:
    /// Attributes ///
    PathHandleGraph& graph;
    IncrementalIdMap<string>& id_map;
    vector<PackedSubgraphOverlay> subgraphs;
    vector<bool> subgraph_partitions;
    unordered_map<nid_t,size_t> node_to_subgraph;
    vector <vector <edge_t> > subgraph_edges;

    HashGraph metagraph;

    /// Methods ///
    Bipartition(PathHandleGraph& graph, IncrementalIdMap<string>& id_map);
    void partition(const unordered_set<nid_t>& node_subset);


};

}

#endif //GFASE_BIPARTITION_HPP
