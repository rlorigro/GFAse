#ifndef GFASE_Hasher22_HPP
#define GFASE_Hasher22_HPP

#include "handlegraph/handle_graph.hpp"
#include "BinarySequence.hpp"
#include "ContactGraph.hpp"
#include "GfaReader.hpp"
#include "Filesystem.hpp"
#include "Sequence.hpp"
#include "CLI11.hpp"
#include "misc.hpp"
#include "spp.h"

using ghc::filesystem::path;
using ghc::filesystem::exists;
using ghc::filesystem::create_directories;
using spp::sparse_hash_set;
using spp::sparse_hash_map;
using handlegraph::HandleGraph;

#include <unordered_set>
#include <map>
#include <iostream>
#include <ostream>
#include <atomic>
#include <thread>
#include <mutex>

using std::unordered_set;
using std::map;
using std::numeric_limits;
using std::stringstream;
using std::ostream;
using std::atomic;
using std::thread;
using std::mutex;
using std::cerr;
using std::min;
using std::max;


namespace gfase{

class HashResult{
public:
    string a;
    string b;
    double ab_over_a;
    double ab_over_b;

    HashResult(const string& a, const string& b, double ab_over_a, double ab_over_b);
    HashResult();
};


// Override the hash function for integer hashing because the hash is provided as the key already
class Hash{
public:
    size_t operator() (uint64_t const& key) const
    {
        return key;
    }
};


class Equal
{
public:
    bool operator() (uint64_t a, uint64_t b) const
    {
        return a == b;
    }
};


// Where to store the names of reads which share hashed-k-mers
using hash_bins_t = vector <unordered_set <string> >;

// Ultimately where the results of LSH are stored
using overlaps_t = sparse_hash_map <string, unordered_map <string, int64_t> >;


class Hasher2{
private:
    // Each bin corresponds to a k-mer, and contains sequence names
    hash_bins_t bins;

    // Can't use a vector for mutexes, and it's probably more efficient to have fewer anyway. So the bin index
    // is downsampled (%) to match the number of mutexes.
    array<mutex,1024> bin_mutexes;

    // The result of counting co-occurring sequences in the hash bins
    overlaps_t overlaps;

    const size_t k;

    size_t n_possible_bins;
    const size_t n_iterations;
    const double total_sample_rate;
    const double iteration_sample_rate;
    const size_t n_threads;

    size_t n_bins;  // Computed from the above values

    // Don't iterate bins with too many hashes
    const size_t max_bin_size = 30;

    // Skip assigning pairs for any match that has fewer than this many hashes
    const size_t min_hashes = 40;

    // How many more bins than the total length of the observed sequence do we want to have, to prevent collisions?
    const size_t bins_scaling_factor = 5;

    static const vector<uint64_t> seeds;

    /// Methods ///
    void hash_sequence(const Sequence& sequence, size_t hash_index);

public:
    Hasher2(size_t k, double sample_rate, size_t n_iterations, size_t n_threads);

    // Main algorithm
    uint64_t hash(const BinarySequence<uint64_t>& kmer, size_t seed_index) const;
    void hash_sequences(const vector<Sequence>& sequences, atomic<size_t>& job_index, const size_t hash_index);
    void hash(const vector<Sequence>& sequences);
    void hash(const HandleGraph& graph, const IncrementalIdMap<string>& id_map);

    // Output/results
    void get_best_matches(map<string, string>& matches, double certainty_threshold) const;
    void get_symmetrical_matches(map<string, string>& symmetrical_matches, double certainty_threshold) const;
    int64_t get_intersection_size(const string& a, const string& b) const;
    void for_each_overlap(const function<void(const pair<string,string>, int64_t weight)>& f) const;
    void for_each_overlap(
            size_t max_hits,
            double min_similarity,
            const function<void(const string& a, const string& b, int64_t n_hashes, int64_t total_hashes)>& f) const;

    // IO
    void write_hash_frequency_distribution() const;
    void write_results(path output_directory) const;

    void convert_to_contact_graph(
            ContactGraph& contact_graph,
            const IncrementalIdMap<string>& id_map,
            double similarity_threshold,
            size_t minimum_hashes,
            size_t max_overlaps) const;

    void deallocate_bins();

};



}

#endif //GFASE_Hasher22_HPP
