#ifndef GFASE_BUBBLEGRAPH_HPP
#define GFASE_BUBBLEGRAPH_HPP

#include "IncrementalIdMap.hpp"
#include "sparsepp/spp.h"
#include "Sam.hpp"

using spp::sparse_hash_map;
using spp::sparse_hash_set;

#include <utility>
#include <vector>
#include <atomic>
#include <mutex>
#include <array>
#include <set>

using std::atomic;
using std::vector;
using std::mutex;
using std::array;
using std::pair;
using std::set;


namespace gfase {

using paired_mappings_t = sparse_hash_map <string, array <set <SamElement>, 2> >;
using unpaired_mappings_t = sparse_hash_map <string, set <SamElement> >;
using contact_map_t = sparse_hash_map <int32_t, sparse_hash_map<int32_t, int32_t> >;


template <class T> class Bubble {
public:
    array<T, 2> ids;
    bool phase;

    Bubble();
    Bubble(T id1, T id2, bool phase);
    void flip();
    T first() const;
    T second() const;
    T get(bool side) const;
    T get_other(T id) const;
    T is_first(T id) const;
    T is_second(T id) const;
    bool contains_id(T id) const;
};


template <class T> Bubble<T>::Bubble(T id1, T id2, bool phase) :
        ids({id1, id2}),
        phase(phase) {}


template <class T> Bubble<T>::Bubble() :
        ids({-1, -1}),
        phase(0) {}


template <class T> void Bubble<T>::flip() {
    phase = not phase;
}


template <class T> T Bubble<T>::get_other(T id) const {
    T other;

    if (ids[0] == id){
        other = ids[1];
    }
    else if (ids[1] == id){
        other = ids[0];
    }
    else{
        throw runtime_error("ERROR: 'get_other': specified id not found in bubble "
                            "(" + to_string(ids[0]) + "," + to_string(ids[1]) + "): " + to_string(id) );
    }

    return other;
}


template <class T> T Bubble<T>::first() const {
    return ids[0 + phase];
}


template <class T> T Bubble<T>::second() const {
    return ids[1 - phase];
}


///
/// s p return
/// 0 1 1
/// 1 0 1
/// 0 0 0
/// 1 1 0
template <class T> T Bubble<T>::get(bool side) const {
    return ids[side != phase];
}


template <class T> T Bubble<T>::is_first(T id) const {
    if (get(0) == id) {
        return true;
    } else {
        return false;
    }
}


template <class T> T Bubble<T>::is_second(T id) const {
    if (get(1) == id) {
        return true;
    } else {
        return false;
    }
}


template <class T> bool Bubble<T>::contains_id(T id) const {
    bool success = false;

    if (get(id) == 0) {
        success = true;
    }

    if (get(id) == 1) {
        success = true;
    }

    return success;
}


// Simple adjacency list that is immutable, generated in one pass from a contact map and shasta bubble naming convention
class BubbleGraph {
private:
    vector <Bubble <int32_t> > bubbles;
    sparse_hash_map <int32_t, int32_t> node_id_to_bubble_id;
    vector <vector <int32_t> > bubble_to_bubble;

    // Additional storage of edges so that they can be iterated quickly/uniquely
    // TODO: understand warning associated with realloc of pair type in sparsepp::sparse_hash_set
    sparse_hash_set <pair <int32_t, int32_t> > bubble_edges;

    void emplace(int32_t id1, int32_t id2, bool phase);

public:
    BubbleGraph();
    BubbleGraph(IncrementalIdMap<string>& id_map, const contact_map_t& contact_map);
    void write_bandage_csv(path output_path, const IncrementalIdMap <string>& id_map) const;
    void generate_bubble_adjacency_from_contact_map(const contact_map_t& contact_map);
    void generate_bubbles_from_shasta_names(IncrementalIdMap <string>& id_map);
    void for_each_adjacent_bubble(int32_t b, const function<void(Bubble<int32_t>& bubble)>& f);
    void for_each_adjacent_bubble(int32_t b, const function<void(const Bubble<int32_t>& bubble)>& f) const;
    void for_each_bubble_edge(const function<void(Bubble<int32_t>& b1, Bubble<int32_t>& b2)>& f);
    void for_each_bubble_edge(const function<void(const Bubble<int32_t>& b1, const Bubble<int32_t>& b2)>& f) const;
    void add_bubble(int32_t node_id_a, int32_t node_id_b);
    int32_t try_add_bubble(int32_t node_id_a, int32_t node_id_b);
    void get_phases(vector<bool>& bubble_phases) const;
    void set_phases(const vector<bool>& bubble_phases);
    void at(size_t i, Bubble<int32_t>& b);
    Bubble<int32_t> at(size_t i) const;
    void for_each_node_id(const function<void(const int32_t id)>& f) const;
    Bubble<int32_t> get_bubble_of_node(int32_t node_id) const;
    int32_t find_bubble_id_of_node(int32_t node_id) const;
    int32_t get_other_side(int32_t node_id) const;
    bool node_is_bubble(int32_t node_id) const;
    void flip(size_t b);
    size_t size() const;
};


void generate_adjacency_matrix(const BubbleGraph& bubbles, const contact_map_t& contact_map, vector <vector <int32_t> >& adjacency);

int64_t compute_total_consistency_score(const BubbleGraph& bubbles, const contact_map_t& contact_map);

int64_t compute_consistency_score(const BubbleGraph& bubbles, size_t bubble_index, const contact_map_t& contact_map);

void random_phase_search(
        const contact_map_t& contact_map,
        const IncrementalIdMap<string>& id_map,
        BubbleGraph bubbles,
        vector<bool>& best_phases,
        atomic<int64_t>& best_score,
        atomic<size_t>& job_index,
        mutex& phase_mutex,
        size_t m_iterations
);

void phase_contacts(
        const contact_map_t& contact_map,
        const IncrementalIdMap<string>& id_map,
        BubbleGraph& bubbles,
        size_t n_threads
);


}


#endif //GFASE_BUBBLEGRAPH_HPP
