#ifndef GFASE_PAIR_HASH_HPP
#define GFASE_PAIR_HASH_HPP

#include <unordered_set>
#include <functional>
#include <utility>

using std::pair;
using std::hash;


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

#endif //GFASE_PAIR_HASH_HPP
