#include "Sam.hpp"
#include "htslib/include/htslib/hts.h"
#include "htslib/include/htslib/sam.h"

#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <istream>
#include <limits>

using std::numeric_limits;
using std::runtime_error;
using std::streamsize;
using std::ifstream;
using std::cerr;
using std::sort;
using std::max;

namespace gfase {


SamElement::SamElement(string& read_name, string& ref_name, uint16_t flag, uint8_t mapq) :
        query_name(read_name),
        ref_name(ref_name),
        flag(flag),
        mapq(mapq)
{}


//SamElement::SamElement(SamElement&& other) noexcept:
//        query_name(std::move(other.query_name)),
//        ref_name(std::move(other.ref_name)),
//        line(other.line),
//        flag(other.flag),
//        mapq(other.mapq)
//{}


//SamElement::SamElement(const SamElement& other) noexcept:
//        query_name(other.query_name),
//        ref_name(other.ref_name),
//        line(other.line),
//        flag(other.flag),
//        mapq(other.mapq)
//{}


SamElement::SamElement() :
        query_name(),
        ref_name(),
        flag(-1),
        mapq(-1)
{}


//SamElement& SamElement::operator=(SamElement&& other) noexcept{
//    if (this != other){
//
//    }
//}


AlignmentBlock::AlignmentBlock(
        int32_t ref_start,
        int32_t ref_stop,
        int32_t query_start,
        int32_t query_stop,
        uint32_t n_matches,
        uint32_t n_mismatches,
        uint32_t n_inserts,
        uint32_t n_deletes,
        bool is_reverse
        ):
    ref_start(ref_start),
    ref_stop(ref_stop),
    query_start(query_start),
    query_stop(query_stop),
    n_matches(n_matches),
    n_mismatches(n_mismatches),
    n_inserts(n_inserts),
    n_deletes(n_deletes),
    is_reverse(is_reverse)
{}


double AlignmentBlock::get_identity() const{
    return double(n_matches) / double(n_matches+n_mismatches+n_inserts+n_deletes+1e-12);
}


char AlignmentBlock::get_reversal_char() const{
    return is_reverse ? '-' : '+';
}


int32_t AlignmentBlock::get_forward_start() const {
    if (is_reverse) {
        return ref_stop;
    } else {
        return ref_start;
    }
}


int32_t AlignmentBlock::get_forward_stop() const {
    if (is_reverse) {
        return ref_start;
    } else {
        return ref_stop;
    }
}


void AlignmentChain::sort_chains(bool by_ref) {
    sort(chain.begin(), chain.end(), [&](const AlignmentBlock& a, const AlignmentBlock& b){
        if (by_ref){
            return a.ref_start < b.ref_start;
        }
        else{
            return a.query_start < b.query_start;
        }
    });
}


void FullAlignmentChain::sort_chains(bool by_ref) {
    sort(chain.begin(), chain.end(), [&](const FullAlignmentBlock& a, const FullAlignmentBlock& b){
        if (by_ref){
            return a.ref_start < b.ref_start;
        }
        else{
            return a.query_start < b.query_start;
        }
    });
}


size_t AlignmentChain::get_approximate_non_overlapping_matches(){
    size_t total_matches = 0;

    auto sorted_chain = *this;
    sorted_chain.sort_chains(true);

    // Don't double count overlaps, use identity to approximate # matches in the overlapping regions
    int32_t prev_stop = -1;
    int32_t prev_start = -1;
    double prev_identity = 0;
    for (auto& c: sorted_chain.chain){
        auto identity = c.get_identity();

        total_matches += c.n_matches;

        if (prev_stop > c.ref_start){
            if (prev_start < c.ref_stop){
                // Skip chains that are entirely contained in previous chains
                // Don't update the prev ref start/stop
                continue;
            }
            else {
                total_matches -= (prev_stop - c.ref_start) * (max(prev_identity, identity));
            }
        }

        prev_stop = c.ref_stop;
        prev_start = c.ref_start;
        prev_identity = identity;
    }

    return total_matches;
}


size_t AlignmentChain::get_total_matches(){
    size_t total_matches = 0;

    for (auto& c: chain){
        total_matches += c.n_matches;
    }

    return total_matches;
}


bool AlignmentChain::empty() {
    return chain.empty();
}


bool FullAlignmentChain::empty() {
    return chain.empty();
}


bool SamElement::is_first_mate() const {
    return (uint16_t(flag) >> 6) & uint16_t(1);
}


bool SamElement::is_second_mate() const {
    return (uint16_t(flag) >> 7) & uint16_t(1);
}


bool SamElement::is_not_primary() const {
    return (uint16_t(flag) >> 8) & uint16_t(1);
}


bool SamElement::is_primary() const {
    return (not is_not_primary());
}


bool SamElement::is_supplementary() const {
    return (uint16_t(flag) >> 11) & uint16_t(1);
}


//// Hash/compare using line number to guarantee unique
//bool sam_comparator(const SamElement& a, const SamElement& b) {
//    return a.line > b.line;
//}


void SamElement::for_each_cigar(const function<void(char type, uint32_t length)>& f) const {
    for (auto& c: cigars){
        char operation = bam_cigar_opchr(c);
        auto length = bam_cigar_oplen(c);
        f(operation, length);
    }
}


void for_element_in_sam_file(path sam_path, const function<void(SamElement& e)>& f){
    ifstream file(sam_path);

    if (not file.is_open() or not file.good()){
        throw runtime_error("ERROR: could not read input file: " + sam_path.string());
    }

    char c;
    size_t n_delimiters = 0;
    int64_t n_lines = 0;

    string read_name;
    string flag_token;
    string ref_name;
    string mapq_token;

    char header_delimiter = '@';

    // If this is a header line, skip to the next line and increment n_lines until we are no longer on a header
    while (file.peek() == header_delimiter){
        cerr << string(1, char(file.peek())) << '\n';
        file.ignore(numeric_limits<streamsize>::max(), '\n');
        n_lines++;
    }

    while (file.get(c)){
        if (c == '\n'){
            n_delimiters = 0;

            SamElement e(read_name, ref_name, int16_t(stoi(flag_token)), int8_t(stoi(mapq_token)));

            f(e);

            read_name.clear();
            flag_token.clear();
            ref_name.clear();
            mapq_token.clear();

            n_lines++;
        }
        else if (c == '\t'){
            n_delimiters++;
        }
        else{
            if (n_delimiters == 0){
                read_name += c;
            }
            else if (n_delimiters == 1){
                flag_token += c;
            }
            else if (n_delimiters == 2){
                ref_name += c;
            }
            else if (n_delimiters == 4){
                mapq_token += c;
            }
        }
    }
}

}


//bool operator<(const gfase::SamElement& a, const gfase::SamElement& b) {
//    return a.line < b.line;
//}
//
//
//bool operator>(const gfase::SamElement& a, const gfase::SamElement& b) {
//    return a.line > b.line;
//}
//
//
//bool operator==(const gfase::SamElement& a, const gfase::SamElement& b) {
//    return a.line == b.line;
//}


ostream& operator<<(ostream& o, const gfase::SamElement& a){
    o << a.query_name << '\n';

    o << '\t' << a.ref_name << '\n';
    o << '\t' << "mapq: " << int(a.mapq) << '\n';
    o << '\t' << "flags: " << bitset<sizeof(a.flag)*8>(a.flag) << " = " << a.flag << '\n';
    o << '\t' << "pair: " << a.is_first_mate() << ' ' << a.is_second_mate() << '\n';
    o << '\t' << "secondary: " << a.is_not_primary() << '\n';
    o << '\t' << "supplementary: " << a.is_supplementary() << '\n';

    return o;
}


ostream& operator<<(ostream& o, const gfase::AlignmentBlock& a){
    o << '\t' << "ref_start: " << a.ref_start << '\n';
    o << '\t' << "ref_stop: " << a.ref_stop << '\n';
    o << '\t' << "query_start: " << a.query_start << '\n';
    o << '\t' << "query_stop: " << a.query_stop << '\n';
    o << '\t' << "n_matches: " << a.n_matches << '\n';

    return o;
}


ostream& operator<<(ostream& o, const gfase::FullAlignmentBlock& a){
    o << '\t' << "ref: " << a.ref_name << '\n';
    o << '\t' << "query: " << a.query_name << '\n';
    o << '\t' << "ref_start: " << a.ref_start << '\n';
    o << '\t' << "ref_stop: " << a.ref_stop << '\n';
    o << '\t' << "query_start: " << a.query_start << '\n';
    o << '\t' << "query_stop: " << a.query_stop << '\n';
    o << '\t' << "n_matches: " << a.n_matches << '\n';
    o << '\t' << "n_mismatches: " << a.n_mismatches << '\n';

    return o;
}

