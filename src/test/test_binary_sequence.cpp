#include "BinarySequence.hpp"
#include "unordered_set"
#include "map"

using gfase::BinarySequence;
using std::unordered_set;
using std::map;
using std::runtime_error;


template <class T> void print_bucket_info(unordered_set <BinarySequence<T> >& set){
    cerr << "total buckets: " << set.bucket_count() << '\n';

    map <size_t,size_t> observed_bucket_sizes;
    for (size_t i=0; i<set.bucket_count(); i++){
        observed_bucket_sizes[set.bucket_size(i)]++;
    }

    for (auto& item: observed_bucket_sizes){
        cerr << item.first << ' ' << item.second << '\n';
    }

}

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


    {
        unordered_set <BinarySequence<uint8_t> > set_u8;
        unordered_set <BinarySequence<int8_t> > set_8;
        unordered_set <BinarySequence<uint16_t> > set_u16;
        unordered_set <BinarySequence<int16_t> > set_16;
        unordered_set <BinarySequence<uint32_t> > set_u32;
        unordered_set <BinarySequence<int32_t> > set_32;
        unordered_set <BinarySequence<uint64_t> > set_u64;
        unordered_set <BinarySequence<int64_t> > set_64;
        unordered_set <BinarySequence<__uint128_t> > set_u128;
        unordered_set <BinarySequence<__int128_t> > set_128;
        unordered_set <BinarySequence<int> > set_int;

        set_u8.reserve(16000);
        set_8.reserve(16000);
        set_u16.reserve(16000);
        set_16.reserve(16000);
        set_u32.reserve(16000);
        set_32.reserve(16000);
        set_u64.reserve(16000);
        set_64.reserve(16000);
        set_u128.reserve(16000);
        set_128.reserve(16000);
        set_int.reserve(16000);

        for (size_t i=0; i<4096; i++) {
            string random_sequence;

            for (size_t j = 0; j < 55; j++) {
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

            set_u8.emplace(bs_u8);
            set_8.emplace(bs_8);
            set_u16.emplace(bs_u16);
            set_16.emplace(bs_16);
            set_u32.emplace(bs_u32);
            set_32.emplace(bs_32);
            set_u64.emplace(bs_u64);
            set_64.emplace(bs_64);
            set_u128.emplace(bs_u128);
            set_128.emplace(bs_128);
            set_int.emplace(bs_int);
        }

        cerr << "set_u8 buckets:" << '\n';
        print_bucket_info(set_u8);
        cerr << "set_8 buckets:" << '\n';
        print_bucket_info(set_8);
        cerr << "set_u16 buckets:" << '\n';
        print_bucket_info(set_u16);
        cerr << "set_16 buckets:" << '\n';
        print_bucket_info(set_16);
        cerr << "set_u32 buckets:" << '\n';
        print_bucket_info(set_u32);
        cerr << "set_32 buckets:" << '\n';
        print_bucket_info(set_32);
        cerr << "set_u64 buckets:" << '\n';
        print_bucket_info(set_u64);
        cerr << "set_64 buckets:" << '\n';
        print_bucket_info(set_64);
        cerr << "set_u128 buckets:" << '\n';
        print_bucket_info(set_u128);
        cerr << "set_128 buckets:" << '\n';
        print_bucket_info(set_128);
        cerr << "set_int buckets:" << '\n';
        print_bucket_info(set_int);
    }

    return 0;
}

