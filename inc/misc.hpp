#ifndef GFASE_MISC_HPP
#define GFASE_MISC_HPP

#include "Filesystem.hpp"
#include "misc.hpp"

using ghc::filesystem::create_directories;
using ghc::filesystem::exists;
using ghc::filesystem::path;

#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <utility>
#include <cstdio>
#include <limits>
#include <vector>
#include <array>
#include <set>
#include <map>


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
using std::getline;
using std::remove;
using std::string;
using std::vector;
using std::array;
using std::pair;
using std::stoi;
using std::pair;
using std::hash;
using std::cerr;
using std::cref;
using std::ref;
using std::set;
using std::map;



namespace gfase{

string join(vector <string> s, char delimiter=' ');

void run_command(string& argument_string);

path align(path output_dir, path ref_path, path query_path, size_t n_threads);

path sam_to_sorted_bam(path sam_path, size_t n_threads, bool remove_sam=true);

void get_query_lengths_from_fasta(path fasta_path, map<string,size_t>& query_lengths);

}


namespace std{

// Code from boost
// Reciprocal of the golden ratio helps spread entropy
//     and handles duplicates.
// See Mike Seymour in magic-numbers-in-boosthash-combine:
//     http://stackoverflow.com/questions/4948780

template <class T>
inline void hash_combine(size_t& seed, T const& v)
{
    seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
}


template <typename A, typename B>
struct hash<pair<A,B> > {
size_t operator()(const pair<A,B>& x) const {
    size_t hash_val = std::hash<A>()(x.first);
    hash_combine(hash_val, x.second);
    return hash_val;
}
};

}

#endif //GFASE_MISC_HPP
