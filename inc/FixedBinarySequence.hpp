#ifndef SIMPLE_DOTPLOT_BINARYSEQUENCE_HPP
#define SIMPLE_DOTPLOT_BINARYSEQUENCE_HPP

#include "MurmurHash3.hpp"
#include "MurmurHash2.hpp"

#include <type_traits>
#include <vector>
#include <bitset>
#include <string>
#include <deque>
#include <array>
#include <iostream>

using std::runtime_error;
using std::is_integral;
using std::to_string;
using std::vector;
using std::bitset;
using std::string;
using std::deque;
using std::array;
using std::cerr;


namespace gfase {


template<class T, size_t T2> class FixedBinarySequence {
public:
    /// Attributes ///
    array <T,T2> sequence;

    static const array<char,4> index_to_base;
    static const array<uint16_t,128> base_to_index;

    /// Methods ///
    FixedBinarySequence();
    FixedBinarySequence(const FixedBinarySequence& s);
    template <class T3> FixedBinarySequence(const T3& s);
    void to_string(string& s, size_t length);
    void print_as_bits() const;
};


template <class T, size_t T2> void FixedBinarySequence<T,T2>::print_as_bits() const {
    for (auto& item: sequence){
        cerr << bitset<sizeof(T)*8>(item) << ' ';
    }
    cerr << '\n';
}


template <class T, size_t T2> bool operator==(const FixedBinarySequence<T,T2>& a, const FixedBinarySequence<T,T2>& b)
{
    return a.sequence == b.sequence;
}


template <class T, size_t T2> const array<char,4> FixedBinarySequence<T,T2>::index_to_base = {'A','C','G','T'};
template <class T, size_t T2> const array<uint16_t,128> FixedBinarySequence<T,T2>::base_to_index = {
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


template<class T, size_t T2> FixedBinarySequence<T,T2>::FixedBinarySequence()
{}


template <class T, size_t T2> FixedBinarySequence<T,T2>::FixedBinarySequence(const FixedBinarySequence<T,T2>& s):
        sequence(s.sequence)
{}


template<class T, size_t T2> template <class T3> FixedBinarySequence<T,T2>::FixedBinarySequence(const T3& s):
    sequence({})  // Bracket initializer ensures the array is filled with zeros
{
    static_assert(is_integral<T>::value, "ERROR: provided type for FixedBinarySequence is not integer");

    size_t length = 0;
    size_t word_index = 0;
    for (auto& c: s){
        T bits = base_to_index.at(c);

        if (bits == 4){
            throw runtime_error("ERROR: non ACGT character encountered in sequence: " + string(c,1));
        }

        uint8_t shift_size = (2*length) % (sizeof(T)*8);

        if (shift_size == 0 and length > 0){
            word_index++;
        }
        else {
            bits <<= shift_size;
        }

//        cerr << bitset<sizeof(T)*8>(sequence.at(word_index)) << '\n';
//        cerr << bitset<sizeof(T)*8>(bits) << '\n';

        sequence.at(word_index) |= bits;
        length++;
    }
}


template<class T, size_t T2> void FixedBinarySequence<T,T2>::to_string(string& s, size_t length){
    if (sequence.empty()){
        return;
    }

    T mask = 3;
    size_t bp_per_word = (sizeof(T)*8)/2;
    T word;

    for (size_t i=0; i<length; i++){
        if (i % bp_per_word == 0) {
            word = sequence[i/bp_per_word];
        }

//        cerr << bitset<sizeof(T)*8>(word) << '\n';
//        cerr << bitset<sizeof(T)*8>(mask) << '\n';

        // Iterate/consume the word and produce bases (chars)
        auto index = word & mask;
        s += index_to_base.at(index);
        word >>= 2;
    }
}


}


namespace std {

template <size_t T2> class hash<gfase::FixedBinarySequence<uint64_t,T2> > {
public:
    size_t operator()(const gfase::FixedBinarySequence<uint64_t,T2>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.sequence.size()*sizeof(uint64_t)), 14741);
    }
};


template <size_t T2> class hash<gfase::FixedBinarySequence<int64_t,T2> > {
public:
    size_t operator()(const gfase::FixedBinarySequence<int64_t,T2>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.sequence.size()*sizeof(int64_t)), 14741);
    }
};


template <size_t T2> class hash<gfase::FixedBinarySequence<uint32_t,T2> > {
public:
    size_t operator()(const gfase::FixedBinarySequence<uint32_t,T2>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.sequence.size()*sizeof(uint32_t)), 14741);
    }
};


template <size_t T2> class hash<gfase::FixedBinarySequence<int32_t,T2> > {
public:
    size_t operator()(const gfase::FixedBinarySequence<int32_t,T2>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.sequence.size()*sizeof(int32_t)), 14741);
    }
};


template <size_t T2> class hash<gfase::FixedBinarySequence<uint16_t,T2> > {
public:
    size_t operator()(const gfase::FixedBinarySequence<uint16_t,T2>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.sequence.size()*sizeof(uint16_t)), 14741);
    }
};


template <size_t T2> class hash<gfase::FixedBinarySequence<int16_t,T2> > {
public:
    size_t operator()(const gfase::FixedBinarySequence<int16_t,T2>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.sequence.size()*sizeof(int16_t)), 14741);
    }
};


template <size_t T2> class hash<gfase::FixedBinarySequence<uint8_t,T2> > {
public:
    size_t operator()(const gfase::FixedBinarySequence<uint8_t,T2>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.sequence.size()*sizeof(uint8_t)), 14741);
    }
};


template <size_t T2> class hash<gfase::FixedBinarySequence<int8_t,T2> > {
public:
    size_t operator()(const gfase::FixedBinarySequence<int8_t,T2>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.sequence.size()*sizeof(int8_t)), 14741);
    }
};


template <size_t T2> class hash<gfase::FixedBinarySequence<__uint128_t,T2> > {
public:
    size_t operator()(const gfase::FixedBinarySequence<__uint128_t,T2>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.sequence.size()*sizeof(__uint128_t)), 14741);
    }
};


template <size_t T2> class hash<gfase::FixedBinarySequence<__int128_t,T2> > {
public:
    size_t operator()(const gfase::FixedBinarySequence<__int128_t,T2>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.sequence.size()*sizeof(__int128_t)), 14741);
    }
};


}



#endif //SIMPLE_DOTPLOT_BINARYSEQUENCE_HPP
