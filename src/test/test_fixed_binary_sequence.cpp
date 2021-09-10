#include "FixedBinarySequence.hpp"
#include <unordered_set>
#include <list>
#include <map>

using gfase::FixedBinarySequence;
using std::unordered_set;
using std::list;
using std::map;
using std::runtime_error;


template <class T, size_t T2> void print_bucket_info(unordered_set <FixedBinarySequence<T,T2> >& set){
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
        FixedBinarySequence<uint64_t,8> bs(s);

        string s2;
        bs.to_string(s2, s.size());

        cerr << s2 << '\n';

        if (s != s2){
            throw runtime_error("ERROR: input sequence does not match output sequence: " + s + " != " + s2);
        }
    }

    {
        FixedBinarySequence<uint16_t,2> bs(s);

        cerr << "Sequence size: " << bs.sequence.size() << '\n';

        string s2;
        bs.to_string(s2, s.size());

        cerr << s2 << '\n';

        s2.clear();
        bs.to_string(s2, s.size());

        cerr << s2 << '\n';

        if (s != s2){
            throw runtime_error("ERROR: input sequence does not match output sequence: " + s + " != " + s2);
        }
    }

    s = "ACCGGGTTTT";
    {
        FixedBinarySequence<uint64_t,2> bs(s);

        string s2;
        bs.to_string(s2, s.size());

        cerr << s2 << '\n';

        if (s != s2){
            throw runtime_error("ERROR: input sequence does not match output sequence: " + s + " != " + s2);
        }
    }

    for (size_t i=0; i<1024; i++){
        string random_sequence;

        for (size_t j=0; j<8; j++){
            size_t base_index = rand() % 4;
            random_sequence += FixedBinarySequence<int,2>::index_to_base[base_index];
        }


        FixedBinarySequence<uint8_t,2> bs_u8(random_sequence);
        string s_bs_u8;
        bs_u8.to_string(s_bs_u8, random_sequence.size());

        if (random_sequence != s_bs_u8){
            throw runtime_error("ERROR: input sequence does not match output sequence:\n\t" + s + "\n\t" + s_bs_u8 + "\nfor type u8");
        }

        FixedBinarySequence<int8_t,2> bs_8(random_sequence);
        string s_bs_8;
        bs_8.to_string(s_bs_8, random_sequence.size());

        if (random_sequence != s_bs_8){
            throw runtime_error("ERROR: input sequence does not match output sequence:\n\t" + s + "\n\t" + s_bs_8 + "\nfor type 8");
        }

        FixedBinarySequence<uint16_t,2> bs_u16(random_sequence);
        string s_bs_u16;
        bs_u16.to_string(s_bs_u16, random_sequence.size());

        if (random_sequence != s_bs_u16){
            throw runtime_error("ERROR: input sequence does not match output sequence:\n\t" + s + "\n\t" + s_bs_u16 + "\nfor type u16");
        }

        FixedBinarySequence<int16_t,2> bs_16(random_sequence);
        string s_bs_16;
        bs_16.to_string(s_bs_16, random_sequence.size());

        if (random_sequence != s_bs_16){
            throw runtime_error("ERROR: input sequence does not match output sequence:\n\t" + s + "\n\t" + s_bs_16 + "\nfor type 16");
        }

        FixedBinarySequence<uint32_t,2> bs_u32(random_sequence);
        string s_bs_u32;
        bs_u32.to_string(s_bs_u32, random_sequence.size());

        if (random_sequence != s_bs_u32){
            throw runtime_error("ERROR: input sequence does not match output sequence:\n\t" + s + "\n\t" + s_bs_u32 + "\nfor type u32");
        }

        FixedBinarySequence<int32_t,2> bs_32(random_sequence);
        string s_bs_32;
        bs_32.to_string(s_bs_32, random_sequence.size());

        if (random_sequence != s_bs_32){
            throw runtime_error("ERROR: input sequence does not match output sequence:\n\t" + s + "\n\t" + s_bs_32 + "\nfor type 32");
        }

        FixedBinarySequence<uint64_t,2> bs_u64(random_sequence);
        string s_bs_u64;
        bs_u64.to_string(s_bs_u64, random_sequence.size());

        if (random_sequence != s_bs_u64){
            throw runtime_error("ERROR: input sequence does not match output sequence:\n\t" + s + "\n\t" + s_bs_u64 + "\nfor type u64");
        }

        FixedBinarySequence<int64_t,2> bs_64(random_sequence);
        string s_bs_64;
        bs_64.to_string(s_bs_64, random_sequence.size());

        if (random_sequence != s_bs_64){
            throw runtime_error("ERROR: input sequence does not match output sequence:\n\t" + s + "\n\t" + s_bs_64 + "\nfor type 64");
        }

        FixedBinarySequence<__uint128_t,2> bs_u128(random_sequence);
        string s_bs_u128;
        bs_u128.to_string(s_bs_u128, random_sequence.size());

        if (random_sequence != s_bs_u128){
            throw runtime_error("ERROR: input sequence does not match output sequence:\n\t" + s + "\n\t" + s_bs_u128 + "\nfor type u128");
        }

        FixedBinarySequence<__int128_t,2> bs_128(random_sequence);
        string s_bs_128;
        bs_128.to_string(s_bs_128, random_sequence.size());

        if (random_sequence != s_bs_128){
            throw runtime_error("ERROR: input sequence does not match output sequence:\n\t" + s + "\n\t" + s_bs_128 + "\nfor type 128");
        }

        FixedBinarySequence<int,2> bs_int(random_sequence);
        string s_bs_int;
        bs_int.to_string(s_bs_int, random_sequence.size());

        if (random_sequence != s_bs_int){
            throw runtime_error("ERROR: input sequence does not match output sequence:\n\t" + s + "\n\t" + s_bs_int + "\nfor type int");
        }
    }

    vector<string> test_kmers = {
            "TTAAAAAAAAAAAAAAATTAAAAAAAAAAAAAAAAAAAAAAAAATAAAAAAAAAA",
            "TTAAAAAAAAAAAAAAATTAAAAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAAA",
            "TTAAACACTTAGCTGAGTTAAACACTTAGCTGAGGCATGGTGATGCATGCCTATA",
            "TTAAATAATTAAAGTCATTAAATAATTAAAGTCATCTTTTCAATGAATGCATTGC",
            "TTAACCCAGTCTCCTTTTTAACCCAGTCTCCTTTGTTAGTTTAGCTAATTTTAGT",
            "TTAATTAAATTTGTACATTAATTAAATTTGTACATCAAAATGATTGAAATAAATC",
            "TTACCAAAATCATAATATTACCAAAATCATAATACATTTAACGTAGACCTGAAAT",
            "TTACTATGAATAATTATTTACTATGAATAATTATATGCCTGCAAATTAGAAAACA",
            "TTACTCCATTCCATATTTTACTCCATTCCATATTACCTCCTATCTCCCACTCTAT",
            "TTAGATGGGTGTGCTAGTTAGATGGGTGTGCTAGCGGGCGCCTGTAATCTCAGCT",
            "TTAGATGGGTGTGCTAGTTAGATGGGTGTGCTAGCGGGTGCCTGTAATCTCAGCT",
            "TTAGCCAAGCATGATGGTTAGCCAAGCATGATGGTGCATGTCTGTGGTCCCAGCT",
            "TTAGCCAGGCATGGTGGTTAGCCAGGCATGGTGGCACTTGCCTGTAATCCCAGCT"
    };

    list <FixedBinarySequence<uint64_t,2> > binary_kmers;
    for (auto& kmer: test_kmers){
        binary_kmers.emplace_back(kmer);
        binary_kmers.back().print_as_bits();
    }


    {
        unordered_set <FixedBinarySequence<uint8_t,16> > set_u8;
        unordered_set <FixedBinarySequence<int8_t,16> > set_8;
        unordered_set <FixedBinarySequence<uint16_t,8> > set_u16;
        unordered_set <FixedBinarySequence<int16_t,8> > set_16;
        unordered_set <FixedBinarySequence<uint32_t,4> > set_u32;
        unordered_set <FixedBinarySequence<int32_t,4> > set_32;
        unordered_set <FixedBinarySequence<uint64_t,2> > set_u64;
        unordered_set <FixedBinarySequence<int64_t,2> > set_64;
        unordered_set <FixedBinarySequence<__uint128_t,1> > set_u128;
        unordered_set <FixedBinarySequence<__int128_t,1> > set_128;
        unordered_set <FixedBinarySequence<int,4> > set_int;

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
                random_sequence += FixedBinarySequence<int,4>::index_to_base[base_index];
            }

            FixedBinarySequence<uint8_t,16> bs_u8(random_sequence);
            FixedBinarySequence<int8_t,16> bs_8(random_sequence);
            FixedBinarySequence<uint16_t,8> bs_u16(random_sequence);
            FixedBinarySequence<int16_t,8> bs_16(random_sequence);
            FixedBinarySequence<uint32_t,4> bs_u32(random_sequence);
            FixedBinarySequence<int32_t,4> bs_32(random_sequence);
            FixedBinarySequence<uint64_t,2> bs_u64(random_sequence);
            FixedBinarySequence<int64_t,2> bs_64(random_sequence);
            FixedBinarySequence<__uint128_t,1> bs_u128(random_sequence);
            FixedBinarySequence<__int128_t,1> bs_128(random_sequence);
            FixedBinarySequence<int,4> bs_int(random_sequence);

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

