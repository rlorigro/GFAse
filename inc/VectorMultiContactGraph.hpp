#ifndef GFASE_VECTORMULTICONTACTGRAPH_HPP
#define GFASE_VECTORMULTICONTACTGRAPH_HPP


#include "IncrementalIdMap.hpp"
#include "MultiContactGraph.hpp"

#include "handlegraph/handle_graph.hpp"
#include "bdsg/hash_graph.hpp"
#include "sparsepp/spp.h"
#include "Filesystem.hpp"
#include "edge.hpp"

using ghc::filesystem::path;
using spp::sparse_hash_map;
using handlegraph::handle_t;
using handlegraph::edge_t;
using handlegraph::nid_t;
using bdsg::HashGraph;


#include <unordered_set>
#include <functional>
#include <utility>
#include <thread>
#include <ostream>
#include <array>
#include <set>
#include <map>
#include <set>

using std::unordered_set;
using std::function;
using std::atomic;
using std::mutex;
using std::pair;
using std::ostream;
using std::array;
using std::map;
using std::set;

namespace gfase {


class VectorMultiNode {
public:
    // To find adjacent nodes
    vector <pair <int32_t,int32_t> > neighbors;

    // To find linked/opposing node in a bubble
    vector<int32_t> alts;

    // Total reads located on this node > mapQ threshold, regardless of pair mapQ
    int64_t coverage;

    // Sequence length of this node
    int32_t length;

    // Which set does this node belong to
    int8_t partition;

    // For simplifying handling of missing nodes
    bool is_null;

    VectorMultiNode();
    VectorMultiNode(int8_t partition);
    VectorMultiNode(const MultiNode& node);
    bool has_alt() const;
};


class VectorMultiContactGraph {
    vector <VectorMultiNode> nodes;
    vector <unordered_map<int32_t,int32_t> > edge_weights;

public:
    // Constructors
    VectorMultiContactGraph()=default;
    VectorMultiContactGraph(const MultiContactGraph& contact_graph);

    // Editing
    void set_partition(int32_t id, int8_t partition);
    void set_partitions(const vector <pair <int32_t,int8_t> >& partitions);

    // Iterating and accessing
    void for_each_edge(const function<void(const pair<int32_t,int32_t>, int32_t weight)>& f) const;
    void get_alt_component(int32_t id, bool validate, alt_component_t& component) const;
    void get_partitions(vector <pair <int32_t,int8_t> >& partitions) const;
    void get_node_ids(vector<int32_t>& ids) const;
    int8_t get_partition(int32_t id) const;
    size_t edge_count(int32_t id) const;
    bool has_alt(int32_t id) const;

    // Optimization
    static double get_score(const VectorMultiNode& a, const VectorMultiNode& b, int32_t weight);
    static double get_score(int8_t p_a, int8_t p_b, int32_t weight);
    double compute_consistency_score(int32_t id) const;
    double compute_consistency_score(int32_t id, int8_t p) const;
    double compute_total_consistency_score() const;
    double compare_total_consistency_score(const MultiContactGraph& other_graph) const;
    void randomize_partitions();
};


}

#endif //GFASE_VECTORMULTICONTACTGRAPH_HPP
