#include "Sam.hpp"

using gfase::AlignmentChain;
using gfase::AlignmentBlock;

#include <iostream>

using std::cerr;


int main(){
    AlignmentChain c;

    size_t stop = 6;
    for (size_t i=0; i < stop - 1; i++){
        AlignmentBlock a;

        a.ref_start = i;
        a.ref_stop = i + 1;

        a.query_start = stop - i - 1;
        a.query_stop = stop - i;

        c.chain.emplace_back(a);
    }

    cerr << "ORIGINAL" << '\n';
    for (auto& item: c.chain){
        cerr << "----" << '\n';
        cerr << item << '\n';
    }
    cerr << '\n';

    cerr << "SORT BY QUERY" << '\n';
    c.sort_chains(false);

    for (auto& item: c.chain){
        cerr << "----" << '\n';
        cerr << item;
    }
    cerr << '\n';

    cerr << "SORT BY REF" << '\n';
    c.sort_chains(true);

    for (auto& item: c.chain){
        cerr << "----" << '\n';
        cerr << item;
    }
    cerr << '\n';

    return 0;
}
