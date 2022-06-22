#ifndef GFASE_FIXEDBINARYSEQUENCE_HPP
#define GFASE_FIXEDBINARYSEQUENCE_HPP

#include "BinarySequence.hpp"
#include "MurmurHash3.hpp"
#include "MurmurHash2.hpp"

#include <type_traits>
#include <ostream>
#include <vector>
#include <bitset>
#include <string>
#include <deque>
#include <array>
#include <cmath>
#include <iostream>

using std::runtime_error;
using std::is_integral;
using std::to_string;
using std::ostream;
using std::vector;
using std::bitset;
using std::string;
using std::deque;
using std::array;
using std::pow;
using std::cerr;


namespace gfase {

/// For any operation involving "length", length must be less than 2*bitcount(T)*T2 because each
/// nucleotide takes two bits to store in binary.
/// \tparam T integer type to use as word storage
/// \tparam T2 number of words to have available for storage in array
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
    template <class T3> FixedBinarySequence(const BinarySequence<T3>& s, size_t fixed_length);
    void get_reverse_complement(FixedBinarySequence<T,T2>& rc, size_t length) const;
    void to_string(string& s, size_t length) const;
    void shift(char c, size_t length);
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


template<class T, size_t T2> FixedBinarySequence<T,T2>::FixedBinarySequence():
        sequence({})
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
            throw runtime_error("ERROR: non ACGT character encountered in sequence: " + string(1,c) + " (ord=" + std::to_string(int(c)) + ")");
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


template<class T, size_t T2> template <class T3> FixedBinarySequence<T,T2>::FixedBinarySequence(const BinarySequence<T3>& s, size_t fixed_length):
    sequence({})  // Bracket initializer ensures the array is filled with zeros
{
    if (sizeof(T) != sizeof(T3)){
        throw runtime_error("ERROR: no implementation for conversion between array of differing word type");
    }
    if (fixed_length < s.length){
        throw runtime_error("ERROR: cannot initialize FixedBinarySequence from longer BinarySequence");
    }

    string test;
    int64_t l = 0;
    for (size_t i=0; i<s.sequence.size(); i++){
        sequence[i] = s.sequence[i];

        l += sizeof(T)*8;

//        cerr << "Converting: " << i << ' ' << l << ' ' << 2*fixed_length << ' ' << (fixed_length*2)%(sizeof(T)*8) << '\n';
//
//        this->to_string(test, fixed_length);
//        cerr << test << '\n';
//        cerr << std::bitset<sizeof(T)*8>(sequence[i]) << '\n';

        if (l > int64_t(fixed_length)*2){
            const T mask = pow(T(2),(fixed_length*2)%(sizeof(T)*8))-1;

//            cerr << std::bitset<sizeof(T)*8>(mask) << '\n';
//            cerr << std::bitset<sizeof(T)*8>(sequence[i]) << '\n';
            sequence[i] &= mask;

//            this->to_string(test, fixed_length);
//            cerr << test << '\n';

            break;
        }
    }
}


/// A function to push a new base but not alter the length of the sequence (as a fixed length queue would)
/// \tparam T
/// \param c
template <class T, size_t T2> void FixedBinarySequence<T,T2>::shift(char c, size_t length){
    T bits = base_to_index.at(c);

    if (bits == 4){
        throw runtime_error("ERROR: non ACGT character encountered in sequence: " + string(1,c) + " (ord=" + std::to_string(int(c)) + ")");
    }

    for (size_t i=0; i<sequence.size(); i++){
        if (i > 0) {
            T mask = 3;
//            cerr << bitset<sizeof(T)*8>(mask) << '\n';
//            cerr << '\n';

            T leftover = sequence[i] & mask;
//            cerr << bitset<sizeof(T)*8>(sequence[i]) << '\n';
//            cerr << bitset<sizeof(T)*8>(leftover) << '\n';

            leftover <<= sizeof(T)*8 - 2;
//            cerr << bitset<sizeof(T)*8>(leftover) << '\n';
//            cerr << '\n';
//
//            cerr << bitset<sizeof(T)*8>(sequence[i-1]) << '\n';
//            cerr << bitset<sizeof(T)*8>(T(pow(2,sizeof(T)*8 - 2) - 1)) << '\n';
            sequence[i-1] &= T(pow(2,sizeof(T)*8 - 2) - 1);
//            cerr << bitset<sizeof(T)*8>(sequence[i-1]) << '\n';

            sequence[i-1] |= leftover;
//            cerr << bitset<sizeof(T)*8>(sequence[i-1]) << '\n';
//
//            cerr << '\n';
//            cerr << '\n';

        }
        sequence[i] >>= 2;
    }

    uint8_t shift_size = (2*(length-1)) % (sizeof(T)*8);

    bits <<= shift_size;

    sequence.back() |= bits;
}


/// Make a new binary sequence with the reverse complement of this one
template <class T, size_t T2> void FixedBinarySequence<T,T2>::get_reverse_complement(FixedBinarySequence<T,T2>& rc, size_t length) const{
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

        // Get the last 2 bits
        auto bits = word & mask;

        // Shift the RC bits over and add the complement of the forward bits
        rc.sequence[(length - i - 1)/bp_per_word] <<= 2;
        rc.sequence[(length - i - 1)/bp_per_word] |= 3 - bits;

        // Advance the bits of the forward complement
        word >>= 2;
    }
}


template<class T, size_t T2> void FixedBinarySequence<T,T2>::to_string(string& s, size_t length) const{
    if (sequence.empty()){
        return;
    }

    s.clear();

    T mask = 3;
    size_t bp_per_word = (sizeof(T)*8)/2;
    T word;

    for (size_t i=0; i<length; i++){
        if (i % bp_per_word == 0) {
            // Move to next word in array
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


template<class T, size_t T2> void get_reverse_complement(const FixedBinarySequence<T,T2>& fc, FixedBinarySequence<T,T2>& rc, size_t length) {
    fc.get_reverse_complement(rc, length);
}


}


namespace std {

template <size_t T2> class hash<gfase::FixedBinarySequence<uint64_t,T2> > {
public:
    size_t operator()(const gfase::FixedBinarySequence<uint64_t,T2>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.sequence.size()*sizeof(uint64_t)), 9223372036854775783);
    }
};


template <size_t T2> class hash<gfase::FixedBinarySequence<int64_t,T2> > {
public:
    size_t operator()(const gfase::FixedBinarySequence<int64_t,T2>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.sequence.size()*sizeof(int64_t)), 9223372036854775783);
    }
};


template <size_t T2> class hash<gfase::FixedBinarySequence<uint32_t,T2> > {
public:
    size_t operator()(const gfase::FixedBinarySequence<uint32_t,T2>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.sequence.size()*sizeof(uint32_t)), 9223372036854775783);
    }
};


template <size_t T2> class hash<gfase::FixedBinarySequence<int32_t,T2> > {
public:
    size_t operator()(const gfase::FixedBinarySequence<int32_t,T2>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.sequence.size()*sizeof(int32_t)), 9223372036854775783);
    }
};


template <size_t T2> class hash<gfase::FixedBinarySequence<uint16_t,T2> > {
public:
    size_t operator()(const gfase::FixedBinarySequence<uint16_t,T2>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.sequence.size()*sizeof(uint16_t)), 9223372036854775783);
    }
};


template <size_t T2> class hash<gfase::FixedBinarySequence<int16_t,T2> > {
public:
    size_t operator()(const gfase::FixedBinarySequence<int16_t,T2>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.sequence.size()*sizeof(int16_t)), 9223372036854775783);
    }
};


template <size_t T2> class hash<gfase::FixedBinarySequence<uint8_t,T2> > {
public:
    size_t operator()(const gfase::FixedBinarySequence<uint8_t,T2>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.sequence.size()*sizeof(uint8_t)), 9223372036854775783);
    }
};


template <size_t T2> class hash<gfase::FixedBinarySequence<int8_t,T2> > {
public:
    size_t operator()(const gfase::FixedBinarySequence<int8_t,T2>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.sequence.size()*sizeof(int8_t)), 9223372036854775783);
    }
};


template <size_t T2> class hash<gfase::FixedBinarySequence<__uint128_t,T2> > {
public:
    size_t operator()(const gfase::FixedBinarySequence<__uint128_t,T2>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.sequence.size()*sizeof(__uint128_t)), 9223372036854775783);
    }
};


template <size_t T2> class hash<gfase::FixedBinarySequence<__int128_t,T2> > {
public:
    size_t operator()(const gfase::FixedBinarySequence<__int128_t,T2>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.sequence.size()*sizeof(__int128_t)), 9223372036854775783);
    }
};


}



#endif //GFASE_FIXEDBINARYSEQUENCE_HPP
