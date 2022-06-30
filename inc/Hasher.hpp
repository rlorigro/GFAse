#ifndef GFASE_HASHER_HPP
#define GFASE_HASHER_HPP

#include "BinarySequence.hpp"
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

#include <unordered_set>
#include <map>
#include <iostream>
#include <ostream>
#include <atomic>
#include <thread>

using std::unordered_set;
using std::map;
using std::numeric_limits;
using std::stringstream;
using std::ostream;
using std::atomic;
using std::thread;
using std::cerr;
using std::min;
using std::max;


namespace gfase{


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
using hash_bins_t = sparse_hash_map <uint64_t, unordered_set <string>, Hash, Equal>;

// How to backtrace for each read to its neighbors
using sketches_t = map <string, sparse_hash_set <uint64_t, Hash, Equal> >;

// Ultimately where the results of LSH are stored
using overlaps_t = sparse_hash_map <string, unordered_map <string, int64_t> >;


class Hasher{
private:
    vector<hash_bins_t> bins_per_iteration;
    vector<sketches_t> sketches_per_iteration;
    overlaps_t overlaps;

    const size_t k;

    size_t n_possible_bins;
    const size_t n_iterations;
    const double total_sample_rate;
    const double iteration_sample_rate;
    const size_t n_threads;

    size_t n_bins;  // Computed from the above values

    static const vector<uint64_t> seeds;

    /// Methods ///
    void hash_sequence(const Sequence& sequence, size_t i);

public:
    Hasher(size_t k, double sample_rate, size_t n_iterations, size_t n_threads);
    static uint64_t hash(const BinarySequence<uint64_t>& kmer, size_t seed_index);
    void write_hash_frequency_distribution(const hash_bins_t& bins) const;
    void hash_sequences(const vector<Sequence>& sequences, atomic<size_t>& job_index);
    void hash(const vector<Sequence>& sequences);
    void get_best_matches(map<string, string>& matches, double certainty_threshold);
    void get_symmetrical_matches(map<string, string>& symmetrical_matches, double certainty_threshold);
    void write_results(path output_directory);
};



}

#endif //GFASE_HASHER_HPP
