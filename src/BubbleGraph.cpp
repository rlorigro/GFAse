#include "BubbleGraph.hpp"
#include "bdsg/internal/hash_map.hpp"
#include "IncrementalIdMap.hpp"
#include "Filesystem.hpp"
#include "sparsepp/spp.h"
#include "misc.hpp"
#include "Sam.hpp"

using gfase::for_element_in_sam_file;
using gfase::IncrementalIdMap;
using gfase::SamElement;

using handlegraph::handle_t;
using handlegraph::nid_t;
using ghc::filesystem::path;
using spp::sparse_hash_map;

#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <stdexcept>
#include <algorithm>
#include <fstream>
#include <utility>
#include <atomic>
#include <thread>
#include <random>
#include <limits>
#include <bitset>
#include <vector>
#include <mutex>
#include <array>
#include <set>

using std::numeric_limits;
using std::unordered_set;
using std::unordered_map;
using std::runtime_error;
using std::streamsize;
using std::exception;
using std::to_string;
using std::function;
using std::ifstream;
using std::ofstream;
using std::shuffle;
using std::string;
using std::vector;
using std::bitset;
using std::thread;
using std::atomic;
using std::array;
using std::mutex;
using std::pair;
using std::stoi;
using std::cerr;
using std::cref;
using std::ref;
using std::set;
using std::min;
using std::max;


namespace gfase {


BubbleGraph::BubbleGraph() :
        bubbles(),
        node_id_to_bubble_id(),
        bubble_to_bubble(),
        bubble_edges() {}


// Strictly adds ids to the id map if pairs (bubble sides) in the shasta convention (.0 or .1 suffix) are incomplete
BubbleGraph::BubbleGraph(IncrementalIdMap<string>& id_map, const contact_map_t& contact_map) :
        bubbles(),
        node_id_to_bubble_id(),
        bubble_to_bubble(),
        bubble_edges() {
    generate_bubbles_from_shasta_names(id_map);
    generate_bubble_adjacency_from_contact_map(contact_map);
}


// Strictly adds ids to the id map if pairs (bubble sides) in the shasta convention (.0 or .1 suffix) are incomplete
BubbleGraph::BubbleGraph(IncrementalIdMap<string>& id_map) :
        bubbles(),
        node_id_to_bubble_id(),
        bubble_to_bubble(),
        bubble_edges() {
    generate_bubbles_from_shasta_names(id_map);
}


// Topology based bubble finding
BubbleGraph::BubbleGraph(const HandleGraph& graph, const contact_map_t& contact_map) :
        bubbles(),
        node_id_to_bubble_id(),
        bubble_to_bubble(),
        bubble_edges() {
    generate_diploid_symmetrical_bubbles_from_graph(graph);
    generate_bubble_adjacency_from_contact_map(contact_map);
}


// Topology based bubble finding
BubbleGraph::BubbleGraph(const HandleGraph& graph) :
        bubbles(),
        node_id_to_bubble_id(),
        bubble_to_bubble(),
        bubble_edges() {
    generate_diploid_symmetrical_bubbles_from_graph(graph);
}


// Load from disk
BubbleGraph::BubbleGraph(path csv_path, IncrementalIdMap<string>& id_map) :
        bubbles(),
        node_id_to_bubble_id(),
        bubble_to_bubble(),
        bubble_edges() {
    generate_bubbles_from_csv(csv_path, id_map);
}


void BubbleGraph::for_each_node_id(const function<void(const int32_t id)>& f) const{
    for (const auto& [node_id, bubble_id]: node_id_to_bubble_id){
        f(node_id);
    }
}


void BubbleGraph::generate_diploid_symmetrical_bubbles_from_graph(const HandleGraph& graph){
    graph.for_each_handle([&](const handle_t& h0) {
        auto id0 = graph.get_id(h0);

        // Do a two-edge walk right/left and left/right
        unordered_set<nid_t> left_second_degree_neighbors;
        unordered_set<nid_t> right_second_degree_neighbors;

        unordered_set<nid_t> left_first_degree_neighbors;
        unordered_set<nid_t> right_first_degree_neighbors;

        graph.follow_edges(h0, true, [&](const handle_t& h1){
            auto id1 = graph.get_id(h1);
            left_first_degree_neighbors.emplace(id1);

            graph.follow_edges(h1, false, [&](const handle_t& h2){
                auto id2 = graph.get_id(h2);

                if (id0 != id2) {
                    left_second_degree_neighbors.emplace(id2);
                }
            });
        });

        graph.follow_edges(h0, false, [&](const handle_t& h1){
            auto id1 = graph.get_id(h1);
            right_first_degree_neighbors.emplace(id1);

            graph.follow_edges(h1, true, [&](const handle_t& h2){
                auto id2 = graph.get_id(h2);

                if (id0 != id2) {
                    right_second_degree_neighbors.emplace(id2);
                }
            });
        });

        int32_t result = -1;

        bool is_symmetrical_bubble = (right_second_degree_neighbors == left_second_degree_neighbors);
        bool is_diploid_bubble = (right_second_degree_neighbors.size() == 1);
        bool is_chainable = (right_first_degree_neighbors.size() < 3 and left_first_degree_neighbors.size() < 3);

        // If there are no second degree neighbors, this unphased subgraph passes
        if (is_symmetrical_bubble and is_diploid_bubble){
            auto id_a = graph.get_id(h0);
            auto id_b = graph.get_id(graph.get_handle(*left_second_degree_neighbors.begin()));

            result = try_add_bubble(int32_t(id_a), int32_t(id_b));
        }
    });
}


void BubbleGraph::generate_bubbles_from_csv(path csv_path, IncrementalIdMap<string>& id_map){
    ifstream file(csv_path);

    if (not (file.is_open() and file.good())){
        throw runtime_error("ERROR: could not write to file: " + csv_path.string());
    }

    string line;
    vector<string> tokens = {""};

    int32_t prev_bubble_id = -1;
    string prev_name;
    bool prev_phase;

    size_t l = 0;
    while (getline(file, line)){
        for (auto& c: line){
            if (c == ','){
                tokens.emplace_back();
            }
            else{
                tokens.back() += c;
            }
        }

        if (l == 0){
            if (tokens[2] != "Phase"){
                throw runtime_error("ERROR: BubbleGraph CSV has incorrect headers");
            }
            tokens = {""};
            l++;
            continue;
        }

        auto& name = tokens[0];
        int32_t bubble_id = stoi(tokens[1]);
        int32_t phase = stoi(tokens[2]);

        if (phase < 0 or phase > 1){
            throw runtime_error("ERROR: BubbleGraph cannot load phase <0 or >1");
        }

        if (bubble_id == prev_bubble_id){
            auto id_a = id_map.try_insert(prev_name);
            auto id_b = id_map.try_insert(name);

            // Phase CSV will indicate phase for one side of the bubble at a time.
            // 0/1 = phase 0
            // 1/0 = phase 1
            if (prev_phase and !phase){
                emplace(int32_t(id_a), int32_t(id_b), 1);
            }
            else if (!prev_phase and phase){
                emplace(int32_t(id_a), int32_t(id_b), 0);
            }
        }

        prev_bubble_id = bubble_id;
        prev_name = name;
        prev_phase = phase;

        tokens = {""};
        l++;
    }
}


// Contact maps are one-sided w.r.t to bubbles. This method simplifies that contact map to a bubble-contact map
void BubbleGraph::generate_bubble_adjacency_from_contact_map(const contact_map_t& contact_map) {
    bubble_to_bubble.resize(bubbles.size());

    for (size_t b = 0; b < bubbles.size(); b++) {
        auto id0 = bubbles[b].get(0);
        auto id1 = bubbles[b].get(1);

        unordered_set <size_t> other_bubbles;

        auto c0 = contact_map.find(id0);
        auto c1 = contact_map.find(id1);

        // Iterate all contacts with id0
        if (c0 != contact_map.end()) {
            for (auto&[other_id, count]: c0->second) {
                auto b_other = find_bubble_id_of_node(other_id);

                // Skip any contacts between regions of the graph that aren't identified as bubbles
                if (b_other < 0){
                    continue;
                }
                other_bubbles.emplace(b_other);
            }
        }

        // Iterate all contacts with id1
        if (c1 != contact_map.end()) {
            for (auto&[other_id, count]: c1->second) {
                auto b_other = find_bubble_id_of_node(other_id);

                // Skip any contacts between regions of the graph that aren't identified as bubbles
                if (b_other < 0){
                    continue;
                }
                other_bubbles.emplace(b_other);
            }
        }

        // Update the adjacency list and the global list of edges
        for (auto& b_other: other_bubbles) {
            bubble_to_bubble[b].emplace_back(b_other);

            // Hash pairs in deterministic ordering
            bubble_edges.emplace(min(b,b_other), max(b,b_other));
        }
    }
}


void BubbleGraph::for_each_adjacent_bubble(int32_t b, const function<void(Bubble<int32_t>& bubble)>& f) {
    for (auto& b_other: bubble_to_bubble[b]) {
        f(bubbles[b_other]);
    }
}


void BubbleGraph::for_each_adjacent_bubble(int32_t b, const function<void(const Bubble<int32_t>& bubble)>& f) const {
    for (const auto& b_other: bubble_to_bubble[b]) {
        f(bubbles[b_other]);
    }
}


void BubbleGraph::for_each_bubble_edge(const function<void(Bubble<int32_t>& b0, Bubble<int32_t>& b1)>& f) {
    for (auto& p: bubble_edges) {
        f(bubbles[p.first], bubbles[p.second]);
    }
}


void BubbleGraph::for_each_bubble_edge(const function<void(const Bubble<int32_t>& b0, const Bubble<int32_t>& b1)>& f) const {
    for (const auto& p: bubble_edges) {
        f(bubbles[p.first], bubbles[p.second]);
    }
}


void BubbleGraph::write_bandage_csv(path output_path, const IncrementalIdMap <string>& id_map) const{
    ofstream file(output_path);

    if (not file.is_open() or not file.good()) {
        throw std::runtime_error("ERROR: could not write to file: " + output_path.string());
    }

    file << "Name" << ',' << "BubbleId" << ',' << "Phase" << ',' << "Color" << '\n';

    for (size_t i=0; i<bubbles.size(); i++) {
        auto b = bubbles.at(i);

        auto id0 = b.get(0);
        auto id1 = b.get(1);

        file << id_map.get_name(id0) << ',' << i << ',' << 0 << ',' << "Dark Orange" << '\n';
        file << id_map.get_name(id1) << ',' << i << ',' << 1 << ',' << "Green Yellow" << '\n';
    }
}


void BubbleGraph::emplace(int32_t id1, int32_t id2, bool phase) {
    node_id_to_bubble_id[id1] = int32_t(bubbles.size());
    node_id_to_bubble_id[id2] = int32_t(bubbles.size());
    bubbles.emplace_back(id1, id2, phase);
}


void BubbleGraph::get_phases(vector<bool>& bubble_phases) const {
    if (bubble_phases.size() != bubbles.size()) {
        bubble_phases.resize(bubbles.size());
    }

    for (size_t b = 0; b < bubbles.size(); b++) {
        bubble_phases[b] = bubbles[b].phase;
    }
}


void BubbleGraph::set_phases(const vector<bool>& bubble_phases) {
    if (bubble_phases.size() != bubbles.size()) {
        throw runtime_error("ERROR: cannot set phase vector with unequal size vector");
    }

    for (size_t b = 0; b < bubbles.size(); b++) {
        bubbles[b].phase = bubble_phases[b];
    }
}


size_t BubbleGraph::size() const {
    return bubbles.size();
}


void BubbleGraph::at(size_t i, Bubble<int32_t>& b) {
    b = bubbles.at(i);
}


Bubble<int32_t> BubbleGraph::at(size_t i) const {
    return bubbles.at(i);
}


Bubble<int32_t> BubbleGraph::get_bubble_of_node(int32_t node_id) const {
    return bubbles.at(node_id_to_bubble_id.at(node_id));
}


int32_t BubbleGraph::find_bubble_id_of_node(int32_t node_id) const {
    auto result = node_id_to_bubble_id.find(node_id);

    if (result != node_id_to_bubble_id.end()){
        return result->second;
    }
    else{
        return -1;
    }
}


bool BubbleGraph::node_is_bubble(int32_t node_id) const {
    return (node_id_to_bubble_id.count(node_id) > 0);
}


int32_t BubbleGraph::get_other_side(int32_t node_id) const{
    auto& bubble = bubbles.at(node_id_to_bubble_id.at(node_id));
    return bubble.get_other(node_id);
}


void BubbleGraph::flip(size_t b) {
    bubbles.at(b).flip();
}


void BubbleGraph::generate_bubbles_from_shasta_names(IncrementalIdMap<string>& id_map) {
    unordered_set <int32_t> visited;

    for (auto&[name, id]: id_map.ids) {
        if (visited.count(int32_t(id)) > 0) {
            continue;
        }

        // Skip any "UR" prefixed nodes
        if (name.empty()){
            continue;
        }
        else if (name[0] == 'U'){
            continue;
        }

        // Split to find last field, which should be 0/1 for shasta PR segments
        auto i = name.find_last_of('.');

        if (i < 2 or i > name.size()) {
            continue;
        }

        auto prefix = name.substr(0, i);
        int64_t side = stoi(name.substr(i + 1));

        // Cheap test to check for proper syntax
        if (side > 1 or side < 0) {
            throw std::runtime_error("ERROR: shasta bubble side not 0/1: " + name);
        }

        // Find complement
        string other_name = prefix + '.' + to_string(1 - side);

        // Look for the other name in the id_map, add it if it doesn't exist
        // Later, will need to assume these missing entries in the contact map are 0
        auto other_id = id_map.try_insert(other_name);
        emplace(int32_t(id), int32_t(other_id), 0);

        visited.emplace(id);
        visited.emplace(other_id);
    }
}


void BubbleGraph::add_bubble(int32_t node_id_a, int32_t node_id_b){
    if (node_id_to_bubble_id.count(node_id_a) != 0){
        throw runtime_error("ERROR: cannot add node id to bubble graph twice: " + to_string(node_id_a));
    }
    if (node_id_to_bubble_id.count(node_id_b) != 0){
        throw runtime_error("ERROR: cannot add node id to bubble graph twice: " + to_string(node_id_b));
    }

    emplace(node_id_a, node_id_b, 0);
}


int32_t BubbleGraph::try_add_bubble(int32_t node_id_a, int32_t node_id_b){
    auto a_result = node_id_to_bubble_id.find(node_id_a);
    auto b_result = node_id_to_bubble_id.find(node_id_b);

    bool a_found = a_result != node_id_to_bubble_id.end();
    bool b_found = b_result != node_id_to_bubble_id.end();

    int32_t bubble_id = -1;

    if (not a_found and not b_found){
        bubble_id = int32_t(bubbles.size());
        emplace(node_id_a, node_id_b, 0);
    }
    else if (a_found and b_found){
        bubble_id = a_result->second;
    }
    else{
        throw runtime_error("ERROR: attempt to create bubble with 1 of 2 IDs that belongs to another bubble: " + (a_found ? to_string(node_id_a) : to_string(node_id_b)));
    }

    return bubble_id;
}


void generate_adjacency_matrix(
        const BubbleGraph& bubbles,
        const contact_map_t& contact_map,
        vector <vector <int32_t> >& adjacency
){

    size_t n = 2*bubbles.size();
    adjacency.resize(n, vector<int32_t>(n, 0));

    for (auto& [id1,map2]: contact_map){
        for (auto& [id2,count]: map2){
            adjacency[id1][id2] = count;
        }
    }
}


int64_t compute_total_consistency_score(
        const BubbleGraph& bubbles,
        const contact_map_t& contact_map
){
    int64_t score = 0;

    bubbles.for_each_bubble_edge([&](const Bubble<int32_t>& b0, const Bubble<int32_t>& b1) {
        auto id0 = b0.get(0);
        auto id1 = b0.get(1);

        auto other_id0 = b1.get(0);
        auto other_id1 = b1.get(1);

        auto r0 = contact_map.find(id0);
        if (r0 != contact_map.end()) {
            auto r00 = r0->second.find(other_id0);
            auto r01 = r0->second.find(other_id1);

            if (r00 != r0->second.end()) {
                score += r00->second;
            }

            if (r01 != r0->second.end()) {
                score -= r01->second;
            }
        }

        auto r1 = contact_map.find(id1);
        if (r1 != contact_map.end()) {
            auto r10 = r1->second.find(other_id0);
            auto r11 = r1->second.find(other_id1);

            if (r10 != r1->second.end()) {
                score -= r10->second;
            }

            if (r11 != r1->second.end()) {
                score += r11->second;
            }
        }
    });

    return score;
}


int64_t compute_consistency_score(
        const BubbleGraph& bubbles,
        size_t bubble_index,
        const contact_map_t& contact_map
){
    int64_t score = 0;

    const Bubble bubble = bubbles.at(int32_t(bubble_index));

    auto id0 = bubble.get(0);
    auto id1 = bubble.get(1);

    bubbles.for_each_adjacent_bubble(int32_t(bubble_index), [&](const Bubble<int32_t>& other_bubble){
        auto other_id0 = other_bubble.get(0);
        auto other_id1 = other_bubble.get(1);

        auto r0 = contact_map.find(id0);
        if (r0 != contact_map.end()){
            auto r00 = r0->second.find(other_id0);
            auto r01 = r0->second.find(other_id1);

            if (r00 != r0->second.end()){
                score += r00->second;
            }

            if (r01 != r0->second.end()){
                score -= r01->second;
            }
        }

        auto r1 = contact_map.find(id1);
        if (r1 != contact_map.end()){
            auto r10 = r1->second.find(other_id0);
            auto r11 = r1->second.find(other_id1);

            if (r10 != r1->second.end()){
                score -= r10->second;
            }

            if (r11 != r1->second.end()){
                score += r11->second;
            }
        }
    });

    return score;
}


void random_phase_search(
        const contact_map_t& contact_map,
        const IncrementalIdMap<string>& id_map,
        BubbleGraph bubbles,
        vector<bool>& best_phases,
        atomic<int64_t>& best_score,
        atomic<size_t>& job_index,
        mutex& phase_mutex,
        size_t m_iterations
){

    size_t m = job_index.fetch_add(1);

    // True random number
    std::random_device rd;

    // Pseudorandom generator with true random seed
    std::mt19937 rng(rd());
    std::uniform_int_distribution<int> uniform_distribution(0,int(bubbles.size()-1));

    int64_t total_score;

    vector<size_t> order(bubbles.size());
    for (size_t i=0; i<order.size(); i++){
        order[i] = i;
    }

    while (m < m_iterations) {
        // Randomly perturb
        for (size_t i=0; i < ((bubbles.size()/10) + 1); i++) {
            bubbles.flip(uniform_distribution(rng));
        }

        for (size_t i=0; i < bubbles.size()*3; i++) {
            auto b = uniform_distribution(rng);

            auto score = compute_consistency_score(bubbles, b, contact_map);
            bubbles.flip(b);
            auto flipped_score = compute_consistency_score(bubbles, b, contact_map);

            // Unflip if original orientation was better
            if (flipped_score < score) {
                bubbles.flip(b);
            }
        }

        total_score = compute_total_consistency_score(bubbles, contact_map);

        phase_mutex.lock();
        if (total_score > best_score) {
            best_score = total_score;
            bubbles.get_phases(best_phases);
        }
        else {
            bubbles.set_phases(best_phases);
        }

        cerr << m << ' ' << best_score << ' ' << total_score << std::flush << '\n';
        phase_mutex.unlock();

        m = job_index.fetch_add(1);
    }
}


void phase_contacts(
        const contact_map_t& contact_map,
        const IncrementalIdMap<string>& id_map,
        BubbleGraph& bubbles,
        size_t n_threads
){

    vector<bool> best_phases(bubbles.size(), false);

    atomic<int64_t> best_score = std::numeric_limits<int64_t>::min();
    atomic<size_t> job_index = 0;

    size_t m_iterations = 10000;

    // Thread-related variables
    vector<thread> threads;
    mutex phase_mutex;

    // Launch threads
    for (uint64_t i=0; i<n_threads; i++){
        try {
            threads.emplace_back(thread(
                    random_phase_search,
                    cref(contact_map),
                    cref(id_map),
                    bubbles,
                    ref(best_phases),
                    ref(best_score),
                    ref(job_index),
                    ref(phase_mutex),
                    m_iterations
            ));
        } catch (const exception &e) {
            cerr << e.what() << "\n";
            exit(1);
        }
    }

    // Wait for threads to finish
    for (auto& t: threads){
        t.join();
    }

    bubbles.set_phases(best_phases);
}

}
