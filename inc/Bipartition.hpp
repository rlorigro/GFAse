#ifndef GFASE_BIPARTITION_HPP
#define GFASE_BIPARTITION_HPP

#include "IncrementalIdMap.hpp"

#include "handlegraph/path_handle_graph.hpp"
#include "handlegraph/handle_graph.hpp"
#include "bdsg/overlays/packed_subgraph_overlay.hpp"
#include "bdsg/hash_graph.hpp"

#include <unordered_set>
#include <functional>
#include <fstream>
#include <array>

using handlegraph::PathHandleGraph;
using handlegraph::HandleGraph;
using handlegraph::handle_t;
using handlegraph::edge_t;
using bdsg::PackedSubgraphOverlay;
using bdsg::HashGraph;

using std::unordered_set;
using std::function;
using std::ostream;
using std::array;


namespace gfase{

class Bipartition {
private:
    /// Attributes ///
    PathHandleGraph& graph;
    IncrementalIdMap<string>& id_map;
    vector<PackedSubgraphOverlay> subgraphs;
    vector<bool> subgraph_partitions;
    unordered_map<nid_t,size_t> node_to_subgraph;
    unordered_map <edge_t, vector<edge_t> > meta_edge_to_edges;
    unordered_set<nid_t> node_subset;

public:
    HashGraph metagraph;

    /// Methods ///
    Bipartition(PathHandleGraph& graph, IncrementalIdMap<string>& id_map, unordered_set<nid_t>& node_subset);
    bool get_partition_of_node(const nid_t& id) const;
    void for_each_subgraph(const function<void(const HandleGraph& subgraph, size_t subgraph_index, bool partition)>& f) const;
    void follow_subgraph_edges(size_t subgraph_index, bool go_left, const function<void(const handle_t& h)>& f);
    void for_each_boundary_node_in_subgraph(size_t subgraph_index, bool left, const function<void(const handle_t& h)>& f);
    void for_each_handle_in_subgraph(size_t subgraph_index, const function<void(const handle_t& h)>& f);
    void partition();
    size_t get_subgraph_size(size_t subgraph_index) const;
    size_t get_subgraph_index(nid_t meta_node) const;
    nid_t get_id(handle_t meta_handle) const;
    void write_parent_graph_csv(ostream& file) const;
    void write_meta_graph_csv(ostream& file) const;
    void print() const;
};

}

#endif //GFASE_BIPARTITION_HPP
