#ifndef GFASE_CONTACTGRAPH_HPP
#define GFASE_CONTACTGRAPH_HPP

#include "IncrementalIdMap.hpp"
#include "BubbleGraph.hpp"          // TODO: Only need `contact_map_t` definition, consider breaking it out

#include "handlegraph/handle_graph.hpp"
#include "bdsg/hash_graph.hpp"
#include "sparsepp/spp.h"
#include "Filesystem.hpp"

using ghc::filesystem::path;
using spp::sparse_hash_map;
using handlegraph::handle_t;
using handlegraph::edge_t;
using handlegraph::nid_t;
using bdsg::HashGraph;


#include <functional>
#include <set>
#include <map>

using std::function;
using std::map;


namespace gfase{


class Node {
public:
    // To find adjacent nodes
    set<int32_t> neighbors;

    // Which set does this node belong to
    int8_t partition;

    Node(int8_t partition);
};


class ContactGraph {
    // Edge map which will only store edges in sorted order {min(a,b), max(a,b)}
    // This is maintained in duplicate with node-level adjacency because there needs to be a 1x fast iteration method
    // to compute total consistency scores
    sparse_hash_map<pair<int32_t,int32_t>, int32_t> edge_weights;
    sparse_hash_map<int32_t,Node> nodes;

public:
    // Constructors
    ContactGraph(const contact_map_t& contact_map, const IncrementalIdMap<string>& id_map);
    ContactGraph()=default;

    // Editing
    void remove_edge(int32_t a, int32_t b);
    void insert_edge(int32_t a, int32_t b);
    void insert_edge(int32_t a, int32_t b, int32_t weight);
    void insert_node(int32_t id, int8_t partition);
    void try_insert_node(int32_t id);
    void insert_node(int32_t id);
    void remove_node(int32_t id);

    // Iterating and accessing
    void for_each_node_neighbor(int32_t node_id, const function<void(int32_t id, const Node& n)>& f);
    void for_each_node(const function<void(int32_t id, const Node& n)>& f);
    void for_each_edge(const function<void(const pair<int32_t,int32_t>, int32_t weight)>& f);
};


}

#endif //GFASE_CONTACTGRAPH_HPP
