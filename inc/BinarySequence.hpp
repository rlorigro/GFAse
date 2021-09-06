#ifndef SIMPLE_DOTPLOT_BINARYSEQUENCE_HPP
#define SIMPLE_DOTPLOT_BINARYSEQUENCE_HPP

#include <type_traits>
#include <vector>
#include <bitset>
#include <string>
#include <array>
#include <iostream>

using std::runtime_error;
using std::is_integral;
using std::vector;
using std::bitset;
using std::string;
using std::array;
using std::cerr;


namespace gfase {


template<class T> class BinarySequence {
public:
    /// Attributes ///
    vector <T> sequence;
    size_t length;

    static const array<char,4> index_to_base;
    static const array<uint16_t,128> base_to_index;

    /// Methods ///
    BinarySequence();
    BinarySequence(string& s);
    void shift(char c);
    void to_string(string& s);
};


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


template<class T> BinarySequence<T>::BinarySequence(string& s):
        length(0)
{
    static_assert(is_integral<T>::value, "ERROR: provided type for BinarySequence is not integer");

    for (auto& c: s){
        shift(c);
    }
}


template <class T> void BinarySequence<T>::shift(char c){
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

//    cerr << int(shift_size) << '\n';
//    cerr << bitset<sizeof(T)*8>(bits) << '\n';
//
//    for (auto& s: sequence) {
//        cerr << bitset<sizeof(T)*8>(s);
//    }
//    cerr << '\n';
}


template<class T> void BinarySequence<T>::to_string(string& s){
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
//            cerr << bitset<sizeof(T)*8>(word) << '\n';

            auto index = word & mask;
            s += index_to_base[index];
            word >>= 2;
        }
    }
}


}


#endif //SIMPLE_DOTPLOT_BINARYSEQUENCE_HPP