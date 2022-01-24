#ifndef GFASE_BINARYSEQUENCE_HPP
#define GFASE_BINARYSEQUENCE_HPP

#include "MurmurHash3.hpp"
#include "MurmurHash2.hpp"

#include <type_traits>
#include <iostream>
#include <vector>
#include <bitset>
#include <string>
#include <deque>
#include <cmath>
#include <array>

using std::runtime_error;
using std::is_integral;
using std::to_string;
using std::vector;
using std::bitset;
using std::string;
using std::deque;
using std::array;
using std::pow;
using std::cerr;


namespace gfase {


template<class T> class BinarySequence {
public:
    /// Attributes ///
    vector <T> sequence;
    uint16_t length;

    static const array<char,4> index_to_base;
    static const array<uint16_t,128> base_to_index;

    /// Methods ///
    BinarySequence();
    BinarySequence(const BinarySequence& s);
    template <class T2> BinarySequence(const T2& s);
    void push_back(char c);
    void shift(char c);
    void to_string(string& s);
    size_t get_byte_length() const;
    void print_as_bits() const;
};


template <class T> void BinarySequence<T>::print_as_bits() const {
    for (auto& item: sequence){
        cerr << bitset<sizeof(T)*8>(item) << ' ';
    }
    cerr << '\n';
}


template <class T> bool operator==(const BinarySequence<T>& a, const BinarySequence<T>& b)
{
    return a.sequence == b.sequence;
}


template <class T> const array<char,4> BinarySequence<T>::index_to_base = {'A','C','G','T'};
template <class T> const array<uint16_t,128> BinarySequence<T>::base_to_index = {
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


template<class T> BinarySequence<T>::BinarySequence():
        length(0)
{}


template <class T> BinarySequence<T>::BinarySequence(const BinarySequence<T>& s):
        sequence(s.sequence),
        length(s.length)
{}


template<class T> template <class T2> BinarySequence<T>::BinarySequence(const T2& s):
        length(0)
{
    static_assert(is_integral<T>::value, "ERROR: provided type for BinarySequence is not integer");

    for (auto& c: s){
        push_back(c);
    }
}


template <class T> void BinarySequence<T>::push_back(char c){
    T bits = base_to_index.at(c);

    if (bits == 4){
        throw runtime_error("ERROR: non ACGT character encountered in sequence: " + string(1,c) + " (ord=" + std::to_string(int(c)) + ")");
    }

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

    if (length >= (1<<sizeof(length)*8) - 1){
        throw runtime_error("ERROR: attempting to append to maximum length BinarySequence: " + std::to_string(length));
    }
}

/// A function to push a new base but not alter the length of the sequence (as a fixed length queue would)
/// \tparam T
/// \param c
template <class T> void BinarySequence<T>::shift(char c){
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


template<class T> void BinarySequence<T>::to_string(string& s){
    if (sequence.empty()){
        return;
    }

    s.clear();

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
//            cerr << bitset<sizeof(T)*8>(word) << '\n';

            auto index = word & mask;
            s += index_to_base[index];
            word >>= 2;
        }
    }
}


template <class T> size_t BinarySequence<T>::get_byte_length() const{
    size_t bit_length = size_t(length)*2;
    size_t byte_length = bit_length/8 + (bit_length % 8 != 0);

    return byte_length;
}


}


namespace std {
template<>
class hash<gfase::BinarySequence<uint64_t> > {
public:
    size_t operator()(const gfase::BinarySequence<uint64_t>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.get_byte_length()), 14741);
    }
};


template<>
class hash<gfase::BinarySequence<int64_t> > {
public:
    size_t operator()(const gfase::BinarySequence<int64_t>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.get_byte_length()), 14741);
    }
};


template<>
class hash<gfase::BinarySequence<uint32_t> > {
public:
    size_t operator()(const gfase::BinarySequence<uint32_t>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.get_byte_length()), 14741);
    }
};


template<>
class hash<gfase::BinarySequence<int32_t> > {
public:
    size_t operator()(const gfase::BinarySequence<int32_t>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.get_byte_length()), 14741);
    }
};


template<>
class hash<gfase::BinarySequence<uint16_t> > {
public:
    size_t operator()(const gfase::BinarySequence<uint16_t>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.get_byte_length()), 14741);
    }
};


template<>
class hash<gfase::BinarySequence<int16_t> > {
public:
    size_t operator()(const gfase::BinarySequence<int16_t>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.get_byte_length()), 14741);
    }
};


template<>
class hash<gfase::BinarySequence<uint8_t> > {
public:
    size_t operator()(const gfase::BinarySequence<uint8_t>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.get_byte_length()), 14741);
    }
};


template<>
class hash<gfase::BinarySequence<int8_t> > {
public:
    size_t operator()(const gfase::BinarySequence<int8_t>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.get_byte_length()), 14741);
    }
};


template<>
class hash<gfase::BinarySequence<__uint128_t> > {
public:
    size_t operator()(const gfase::BinarySequence<__uint128_t>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.get_byte_length()), 14741);
    }
};


template<>
class hash<gfase::BinarySequence<__int128_t> > {
public:
    size_t operator()(const gfase::BinarySequence<__int128_t>& s) const {
        return MurmurHash64A(s.sequence.data(), int(s.get_byte_length()), 14741);
    }
};


}



#endif //GFASE_BINARYSEQUENCE_HPP
