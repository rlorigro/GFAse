#ifndef GFASE_BIPARTITION_HPP
#define GFASE_BIPARTITION_HPP

#include "IncrementalIdMap.hpp"

#include "sparsepp/spp.h"

#include "handlegraph/path_handle_graph.hpp"
#include "handlegraph/handle_graph.hpp"
#include "bdsg/overlays/packed_subgraph_overlay.hpp"
#include "bdsg/hash_graph.hpp"

#include <unordered_set>
#include <functional>
#include <fstream>
#include <array>

using bdsg::PathHandleGraph;
using bdsg::HandleGraph;
using bdsg::handle_t;
using bdsg::edge_t;
using bdsg::PackedSubgraphOverlay;
using bdsg::HashGraph;
using spp::sparse_hash_map;

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
    sparse_hash_map<size_t,PackedSubgraphOverlay> subgraphs;
    sparse_hash_map<size_t,bool> subgraph_partitions;
    sparse_hash_map<nid_t,size_t> node_to_subgraph;
    sparse_hash_map <edge_t, unordered_set<edge_t> > meta_edge_to_edges;
    unordered_set<nid_t> node_subset;
    size_t max_id = 0;

public:
    HashGraph metagraph;

    /// Methods ///
    Bipartition(PathHandleGraph& graph, IncrementalIdMap<string>& id_map, unordered_set<nid_t>& node_subset);
    bool get_partition_of_node(const nid_t& id) const;
    bool get_partition_of_subgraph(const size_t subgraph_index) const;
    void for_each_subgraph(const function<void(const HandleGraph& subgraph, size_t subgraph_index, bool partition)>& f) const;
    void follow_subgraph_edges(size_t subgraph_index, bool go_left, const function<void(const handle_t& h)>& f);
    void for_each_boundary_node_in_subgraph(size_t subgraph_index, bool left, const function<void(const handle_t& h)>& f);
    void for_each_handle_in_subgraph(size_t subgraph_index, const function<void(const handle_t& h)>& f);
    nid_t get_id_of_parent_handle(const handle_t& h);
    nid_t get_id_of_parent_handle(const string& name);
    string get_name_of_parent_node(nid_t id);
    void merge_subgraphs(size_t subgraph_index_a, size_t subgraph_index_b);
    void partition();
    size_t get_subgraph_size(size_t subgraph_index) const;
    size_t get_subgraph_index_of_parent_node(nid_t id) const;
    nid_t get_id(handle_t meta_handle) const;
    void write_parent_graph_csv(ostream& file) const;
    void write_meta_graph_csv(ostream& file) const;
    void print() const;
    size_t size() const;
};

}

#endif //GFASE_BIPARTITION_HPP
