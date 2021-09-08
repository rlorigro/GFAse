#ifndef SIMPLE_DOTPLOT_BINARYSEQUENCE_HPP
#define SIMPLE_DOTPLOT_BINARYSEQUENCE_HPP

#include "MurmurHash3.hpp"
#include "MurmurHash2.hpp"

#include <type_traits>
#include <vector>
#include <bitset>
#include <string>
#include <array>
#include <iostream>

using std::runtime_error;
using std::is_integral;
using std::to_string;
using std::vector;
using std::bitset;
using std::string;
using std::array;
using std::cerr;


namespace gfase {


template<class T> class FixedBinarySequence {
public:
    /// Attributes ///
    vector <T> sequence;

    static const array<char,4> index_to_base;
    static const array<uint16_t,128> base_to_index;

    /// Methods ///
    FixedBinarySequence();
    FixedBinarySequence(string& s);
    void to_string(string& s, size_t length);
    void print_as_bits() const;
};


template <class T> void FixedBinarySequence<T>::print_as_bits() const {
    for (auto& item: sequence){
        cerr << bitset<sizeof(T)*8>(item) << ' ';
    }
    cerr << '\n';
}


template <class T> bool operator==(const FixedBinarySequence<T>& a, const FixedBinarySequence<T>& b)
{
    return a.sequence == b.sequence;
}


template <class T> const array<char,4> FixedBinarySequence<T>::index_to_base = {'A','C','G','T'};
template <class T> const array<uint16_t,128> FixedBinarySequence<T>::base_to_index = {
        4,4,4,4,4,4,4,4,4,4,      // 0
        4,4,4,4,4,4,4,4,4,4,      // 10
        4,4,4,4,4,4,4,4,4,4,      // 20
        4,4,4,4,4,4,4,4,4,4,      // 30
        4,4,4,4,4,4,4,4,4,4,      // 40
        4,4,4,4,4,4,4,4,4,4,      // 50
        4,4,4,4,4,0,4,1,4,4,      // 60  A = 65, C = 67
        4,2,4,4,4,4,4,4,4,4,      // 70  G = 71
        4,4,4,4,3,4,4,4,4,4,      // 80  T = 84
        4,4,4,4,4,4,4,4,4,4,      // 90
        4,4,4,4,4,4,4,4,4,4,      // 100
        4,4,4,4,4,4,4,4,4,4,      // 110
        4,4,4,4,4,4,4,4           // 120
};


template<class T> FixedBinarySequence<T>::FixedBinarySequence()
{}


template<class T> FixedBinarySequence<T>::FixedBinarySequence(string& s)
{
    static_assert(is_integral<T>::value, "ERROR: provided type for FixedBinarySequence is not integer");

    size_t length = 0;
    for (auto& c: s){
        T bits = base_to_index.at(c);
        uint8_t shift_size = (2*length) % (sizeof(T)*8);

        // If we have reached the beginning of a new word, append the word vector with 0
        if (shift_size == 0){
            sequence.emplace_back(0);
        }
        // Otherwise shift the last word over to prepare to insert another base
        else {
            bits <<= shift_size;
        }

        sequence.back() |= bits;
        length++;
    }
}


template<class T> void FixedBinarySequence<T>::to_string(string& s, size_t length){
    if (sequence.empty()){
        return;
    }

    T mask = 3;

    size_t n_bits = sizeof(T)*8;
    size_t n_leftover_bases = length - (n_bits*(sequence.size() - 1))/2;

    for (size_t i=0; i<sequence.size(); i++){
        T word = sequence[i];
        size_t l = n_bits/2;

        if (i == sequence.size() - 1){
            l = n_leftover_bases;
        }

        // Iterate/consume the word and produce bases (chars)
        for (size_t j=0; j<l; j++){
            auto index = word & mask;
            s += index_to_base[index];
            word >>= 2;
        }
    }
}


}


namespace std {
template<>
class hash<gfase::FixedBinarySequence<uint64_t> > {
public:
    size_t operator()(const gfase::FixedBinarySequence<uint64_t>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.sequence.size()*sizeof(uint64_t)), 14741);
    }
};


template<>
class hash<gfase::FixedBinarySequence<int64_t> > {
public:
    size_t operator()(const gfase::FixedBinarySequence<int64_t>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.sequence.size()*sizeof(int64_t)), 14741);
    }
};


template<>
class hash<gfase::FixedBinarySequence<uint32_t> > {
public:
    size_t operator()(const gfase::FixedBinarySequence<uint32_t>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.sequence.size()*sizeof(uint32_t)), 14741);
    }
};


template<>
class hash<gfase::FixedBinarySequence<int32_t> > {
public:
    size_t operator()(const gfase::FixedBinarySequence<int32_t>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.sequence.size()*sizeof(int32_t)), 14741);
    }
};


template<>
class hash<gfase::FixedBinarySequence<uint16_t> > {
public:
    size_t operator()(const gfase::FixedBinarySequence<uint16_t>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.sequence.size()*sizeof(uint16_t)), 14741);
    }
};


template<>
class hash<gfase::FixedBinarySequence<int16_t> > {
public:
    size_t operator()(const gfase::FixedBinarySequence<int16_t>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.sequence.size()*sizeof(int16_t)), 14741);
    }
};


template<>
class hash<gfase::FixedBinarySequence<uint8_t> > {
public:
    size_t operator()(const gfase::FixedBinarySequence<uint8_t>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.sequence.size()*sizeof(uint8_t)), 14741);
    }
};


template<>
class hash<gfase::FixedBinarySequence<int8_t> > {
public:
    size_t operator()(const gfase::FixedBinarySequence<int8_t>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.sequence.size()*sizeof(int8_t)), 14741);
    }
};


template<>
class hash<gfase::FixedBinarySequence<__uint128_t> > {
public:
    size_t operator()(const gfase::FixedBinarySequence<__uint128_t>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.sequence.size()*sizeof(__uint128_t)), 14741);
    }
};


template<>
class hash<gfase::FixedBinarySequence<__int128_t> > {
public:
    size_t operator()(const gfase::FixedBinarySequence<__int128_t>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.sequence.size()*sizeof(__int128_t)), 14741);
    }
};


}



#endif //SIMPLE_DOTPLOT_BINARYSEQUENCE_HPP
