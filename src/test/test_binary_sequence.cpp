#include "BinarySequence.hpp"

using gfase::BinarySequence;
using std::runtime_error;


int main(){
    string s = "ATATATATATAT";

    {
        BinarySequence<uint64_t> bs(s);

        string s2;
        bs.to_string(s2);

        cerr << s2 << '\n';

        if (s != s2){
            throw runtime_error("ERROR: input sequence does not match output sequence: " + s + " != " + s2);
        }
    }

    {
        BinarySequence<uint16_t> bs(s);

        cerr << "Sequence size: " << bs.sequence.size() << '\n';

        string s2;
        bs.to_string(s2);

        cerr << s2 << '\n';

        s2.clear();
        bs.to_string(s2);

        cerr << s2 << '\n';

        if (s != s2){
            throw runtime_error("ERROR: input sequence does not match output sequence: " + s + " != " + s2);
        }
    }

    s = "ACCGGGTTTT";
    {
        BinarySequence<uint64_t> bs(s);

        string s2;
        bs.to_string(s2);

        cerr << s2 << '\n';

        if (s != s2){
            throw runtime_error("ERROR: input sequence does not match output sequence: " + s + " != " + s2);
        }
    }

    for (size_t i=0; i<1024; i++){
        string random_sequence;

        for (size_t j=0; j<i; j++){
            size_t base_index = rand() % 4;
            random_sequence += BinarySequence<int>::index_to_base[base_index];
        }

        BinarySequence<uint8_t> bs_u8(random_sequence);
        BinarySequence<int8_t> bs_8(random_sequence);
        BinarySequence<uint16_t> bs_u16(random_sequence);
        BinarySequence<int16_t> bs_16(random_sequence);
        BinarySequence<uint32_t> bs_u32(random_sequence);
        BinarySequence<int32_t> bs_32(random_sequence);
        BinarySequence<uint64_t> bs_u64(random_sequence);
        BinarySequence<int64_t> bs_64(random_sequence);
        BinarySequence<__uint128_t> bs_u128(random_sequence);
        BinarySequence<__int128_t> bs_128(random_sequence);
        BinarySequence<int> bs_int(random_sequence);

        string s_bs_u8;
        bs_u8.to_string(s_bs_u8);

        if (random_sequence != s_bs_u8){
            throw runtime_error("ERROR: input sequence does not match output sequence:\n\t" + s + "\n\t" + s_bs_u8 + "\nfor type u8");
        }

        string s_bs_8;
        bs_8.to_string(s_bs_8);

        if (random_sequence != s_bs_8){
            throw runtime_error("ERROR: input sequence does not match output sequence:\n\t" + s + "\n\t" + s_bs_8 + "\nfor type 8");
        }

        string s_bs_u16;
        bs_u16.to_string(s_bs_u16);

        if (random_sequence != s_bs_u16){
            throw runtime_error("ERROR: input sequence does not match output sequence:\n\t" + s + "\n\t" + s_bs_u16 + "\nfor type u16");
        }

        string s_bs_16;
        bs_16.to_string(s_bs_16);

        if (random_sequence != s_bs_16){
            throw runtime_error("ERROR: input sequence does not match output sequence:\n\t" + s + "\n\t" + s_bs_16 + "\nfor type 16");
        }

        string s_bs_u32;
        bs_u32.to_string(s_bs_u32);

        if (random_sequence != s_bs_u32){
            throw runtime_error("ERROR: input sequence does not match output sequence:\n\t" + s + "\n\t" + s_bs_u32 + "\nfor type u32");
        }

        string s_bs_32;
        bs_32.to_string(s_bs_32);

        if (random_sequence != s_bs_32){
            throw runtime_error("ERROR: input sequence does not match output sequence:\n\t" + s + "\n\t" + s_bs_32 + "\nfor type 32");
        }

        string s_bs_u64;
        bs_u64.to_string(s_bs_u64);

        if (random_sequence != s_bs_u64){
            throw runtime_error("ERROR: input sequence does not match output sequence:\n\t" + s + "\n\t" + s_bs_u64 + "\nfor type u64");
        }

        string s_bs_64;
        bs_64.to_string(s_bs_64);

        if (random_sequence != s_bs_64){
            throw runtime_error("ERROR: input sequence does not match output sequence:\n\t" + s + "\n\t" + s_bs_64 + "\nfor type 64");
        }

        string s_bs_u128;
        bs_u128.to_string(s_bs_u128);

        if (random_sequence != s_bs_u128){
            throw runtime_error("ERROR: input sequence does not match output sequence:\n\t" + s + "\n\t" + s_bs_u128 + "\nfor type u128");
        }

        string s_bs_128;
        bs_128.to_string(s_bs_128);

        if (random_sequence != s_bs_128){
            throw runtime_error("ERROR: input sequence does not match output sequence:\n\t" + s + "\n\t" + s_bs_128 + "\nfor type 128");
        }

        string s_bs_int;
        bs_int.to_string(s_bs_int);

        if (random_sequence != s_bs_int){
            throw runtime_error("ERROR: input sequence does not match output sequence:\n\t" + s + "\n\t" + s_bs_int + "\nfor type int");
        }
    }

//    {
//        BinarySequence<float> bs(s);
//
//        string s2;
//        bs.to_string(s2);
//
//        cerr << s2 << '\n';
//
//        if (s != s2){
//            throw runtime_error("ERROR: input sequence does not match output sequence: " + s + " != " + s2);
//        }
//    }


    return 0;
}

