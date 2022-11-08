#ifndef GFASE_MULTICONTACTGRAPH_HPP
#define GFASE_MULTICONTACTGRAPH_HPP

#include "IncrementalIdMap.hpp"
#include "BubbleGraph.hpp"          // TODO: Only need `contact_map_t` definition, consider breaking it out

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


namespace gfase{


using alt_component_t = pair <set<int32_t>, set<int32_t> >;


class NonBipartiteEdgeException : public std::runtime_error{
public:
    alt_component_t component_a;
    alt_component_t component_b;
    int32_t a;
    int32_t b;
    vector<int32_t> conflicts_0;
    vector<int32_t> conflicts_1;
    string message;

    [[nodiscard]] const char* what() const noexcept override;

    NonBipartiteEdgeException(const alt_component_t& c_a, const alt_component_t& c_b, int32_t a, int32_t b);
};


class MultiNode {
public:
    // To find adjacent nodes
    set<int32_t> neighbors;

    // To find linked/opposing node in a bubble
    set<int32_t> alts;

    // Total reads located on this node > mapQ threshold, regardless of pair mapQ
    int64_t coverage;

    // Sequence length of this node
    int32_t length;

    // Which set does this node belong to
    int8_t partition;

    MultiNode(int8_t partition);
    bool has_alt() const;
};


class MultiContactGraph {
    // Edge map which will only store edges in sorted order {min(a,b), max(a,b)}
    // This is maintained in duplicate with node-level adjacency because there needs to be a 1x fast iteration method
    // to compute total consistency scores
    unordered_map<pair<int32_t,int32_t>, int32_t> edge_weights;
    unordered_map<int32_t,MultiNode> nodes;

    static const array<string,3> colors;
    int32_t max_id;

    // No safety checks built in, only should be called when it's known that the nodes exist and the edge does not.
    void insert_edge(int32_t a, int32_t b, int32_t weight);

public:
    // Constructors
    MultiContactGraph(const contact_map_t& contact_map, const IncrementalIdMap<string>& id_map);
    MultiContactGraph(path csv_path, const IncrementalIdMap<string>& id_map);
    MultiContactGraph();

    // Editing
    void remove_edge(int32_t a, int32_t b);
    void try_insert_edge(int32_t a, int32_t b);
    void try_insert_edge(int32_t a, int32_t b, int32_t weight);
    void increment_edge_weight(int32_t a, int32_t b, int32_t value);
    void increment_coverage(int32_t a, int64_t value);
    void set_node_coverage(int32_t a, int64_t value);
    void set_node_length(int32_t a, int32_t length);
    void insert_node(int32_t id);
    void insert_node(int32_t id, int8_t partition);
    void try_insert_node(int32_t id);
    void try_insert_node(int32_t id, int8_t partition);
    void remove_node(int32_t id);
    void set_partition(int32_t id, int8_t partition);
    void set_partition(const alt_component_t& component, int8_t partition);
    void add_alt(int32_t a, int32_t b);
    void add_alt(const alt_component_t& a, const alt_component_t& b, bool remove_weights=true);
    size_t edge_count(int32_t id) const;
    size_t edge_count() const;

    // Iterating and accessing
    void for_each_node_neighbor(int32_t id, const function<void(int32_t id_other, const MultiNode& n)>& f) const;
    void for_each_node_neighbor(int32_t id, const function<void(int32_t id_other)>& f) const;
    void for_each_node(const function<void(int32_t id, const MultiNode& n)>& f) const;
    void for_each_node(const function<void(int32_t id)>& f) const;
    void for_each_edge(const function<void(const pair<int32_t,int32_t>, int32_t weight)>& f) const;
    void for_each_alt(int32_t id, const function<void(int32_t alt_id, MultiNode& alt_node)>& f);
    void for_each_double_alt(int32_t id, const function<void(int32_t alt_id, MultiNode& alt_node)>& f);
    void get_node_ids(vector<int32_t>& ids);
    bool has_alt(int32_t id) const;
    bool has_node(int32_t id) const;
    bool has_edge(int32_t a, int32_t b) const;
    int64_t get_node_coverage(int32_t id) const;
    int32_t get_node_length(int32_t id) const;
    int32_t get_edge_weight(int32_t a, int32_t b) const;
    void get_alt_component(int32_t id, bool validate, alt_component_t& component) const;
    void get_alt_components(vector <alt_component_t>& alt_components) const;
    void get_alt_component_representatives(vector<int32_t>& representative_ids) const;
    int8_t get_partition(int32_t id) const;

    // Optimization
    double get_score(const MultiNode& a, const MultiNode& b, int32_t weight) const;
    double get_score(int32_t id_a, int32_t id_b) const;
    double compute_total_consistency_score() const;
    double compute_consistency_score(int32_t id) const;
    double compute_consistency_score(alt_component_t& component) const;
    void get_partitions(vector <pair <int32_t,int8_t> >& partitions) const;
    void set_partitions(const vector <pair <int32_t,int8_t> >& partitions);
    void randomize_partitions();

    // IO
    void write_contact_map(path output_path, const IncrementalIdMap<string>& id_map) const;
    void write_bandage_csv(path output_path, const IncrementalIdMap<string>& id_map) const;
    void write_node_data(path output_path, const IncrementalIdMap<string>& id_map) const;

    // Misc
    void assert_component_is_valid(const alt_component_t& component) const;
    bool components_are_compatible(const alt_component_t& component_a, const alt_component_t& component_b) const;
    bool of_same_component(int32_t id_a, int32_t id_b) const;
    bool of_same_component_side(int32_t id_a, int32_t id_b) const;
    void validate_alts();
    size_t size() const;
    size_t get_max_id() const;

    void merge_components(
        const alt_component_t& component_a,
        const alt_component_t& component_b,
        alt_component_t& merged_component) const;

    void get_alts_from_shasta_names(const IncrementalIdMap<string>& id_map);

};


ostream& operator<<(ostream& o, const gfase::MultiNode& n);

}

#endif //GFASE_MULTICONTACTGRAPH_HPP
