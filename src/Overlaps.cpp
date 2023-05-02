#include "Overlaps.hpp"

#include <stdlib.h>
#include <stdexcept>

using std::to_string;
using std::make_pair;
using std::runtime_error;
using std::move;

///
/// M   I   D   N   S   H   P   =   X
/// 0   1   2   3   4   5   6   7   8
///


    // 0   1   2   3   4   5   6   7   8   9
const array<uint8_t, 128> CigarOperation::cigar_code = {
    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,  // 0
    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,  // 10
    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,  // 20
    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,  // 30
    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,  // 40
    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,  // 50
    10, 7,  10, 10, 10, 10, 10, 10, 2,  10,  // 60    =61, D68
    10, 10, 5,  1,  10, 10, 10, 0,  3,  10,  // 70    H72, I73, M77, N78
    6,  10, 10, 4,  10, 10, 10, 10, 8,  10,  // 80    P80, S83, X88
    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,  // 90
    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,  // 100
    10, 10, 10, 10, 10, 10, 10, 10, 10, 10,  // 110
    10, 10, 10, 10, 10, 10, 10, 10};         // 120


const array<bool, 9> CigarOperation::is_ref_move = {
    true,      //MATCH      0
    false,     //INS        1
    true,      //DEL        2
    true,      //REF_SKIP   3
    false,     //SOFT_CLIP  4
    false,     //HARD_CLIP  5
    false,     //PAD        6
    true,      //EQUAL      7
    true};     //MISMATCH   8


const array<bool, 9> CigarOperation::is_query_move = {
    true,     //MATCH      0
    true,     //INS        1
    false,    //DEL        2
    false,    //REF_SKIP   3
    true,     //SOFT_CLIP  4
    false,    //HARD_CLIP  5
    false,    //PAD        6
    true,     //EQUAL      7
    true};    //MISMATCH   8

const array<char, 9> CigarOperation::cigar_type = {'M','I','D','N','S','H','P','=','X'};


CigarOperation::CigarOperation(uint32_t length, char type):
    code(cigar_code[type]),
    op_length(length)
{
    if (code > 8){
        throw runtime_error("ERROR: unrecognized cigar character: " + string(1,type) + " has ASCII value: " + to_string(int(code)));
    }
}

char CigarOperation::type() const {
    return cigar_type[code];
}

uint32_t CigarOperation::length() const {
    return op_length;
}

uint32_t CigarOperation::ref_length() const {
    return is_ref_move[code] ? op_length : 0;
}

uint32_t CigarOperation::query_length() const {
    return is_query_move[code] ? op_length : 0;
}

Cigar::Cigar(const std::string& cigar_string) {
    
    if (cigar_string == "*") {
        // we handle the placeholder for an empy cigar from GFA
        return;
    }
    
    string token;
    for (auto& c: cigar_string) {
        if (isdigit(c)) {
            token += c;
        }
        else {
            auto len = stol(token);
            if (len > 0) {
                operations.emplace_back(len, c);
                token.clear();
            }
        }
    }
}

string Cigar::get_string() const {
    
    string s;
    
    // Count up the cigar operations
    for (const auto& c: operations) {
        s += to_string(c.length());
        s += c.type();
    }
    
    if (s.empty()) {
        s = "0M";
    }
    
    return s;
}

Cigar Cigar::reverse() const {
    Cigar reversed;
    for (auto it = operations.rbegin(); it != operations.rend(); ++it) {
        if (it->type() == 'S' || it->type() == 'H' || it->type() == 'N' || it->type() == 'P') {
            throw runtime_error("Clipped or unaligned CIGARs cannot be trivially reversed");
        }
        else if (it->type() == 'I') {
            reversed.operations.emplace_back(it->length(), 'D');
        }
        else if (it->type() == 'D') {
            reversed.operations.emplace_back(it->length(), 'I');
        }
        else {
            reversed.operations.emplace_back(it->length(), it->type());
        }
    }
    return reversed;
}

pair<size_t, size_t> Cigar::aligned_length() const {
    pair<size_t, size_t> lengths(0, 0);
    for (const auto& op : operations) {
        lengths.first += op.ref_length();
        lengths.second += op.query_length();
    }
    return lengths;
}

bool Cigar::empty() const {
    return operations.empty();
}

size_t Cigar::size() const {
    return operations.size();
}

const CigarOperation& Cigar::at(size_t i) const {
    return operations.at(i);
}

Cigar::iterator Cigar::begin() const {
    return operations.begin();
}

Cigar::iterator Cigar::end() const {
    return operations.end();
}

void Overlaps::record_overlap(const HandleGraph& graph, handle_t a, handle_t b, const string& cigar) {
    record_overlap(graph, a, b, Cigar(cigar));
}

void Overlaps::record_overlap(const HandleGraph& graph, handle_t a, handle_t b, const Cigar& cigar) {
    if (!cigar.empty()) {
        overlaps[graph.edge_handle(a, b)] = cigar;
    }
}

void Overlaps::remove_overlap(const HandleGraph& graph, handle_t a, handle_t b) {
    auto it = overlaps.find(graph.edge_handle(a, b));
    if (it != overlaps.end()) {
        overlaps.erase(it);
    }
}


bool Overlaps::has_overlap(const HandleGraph& graph, handle_t a, handle_t b) const {
    return overlaps.count(graph.edge_handle(a, b));
}

Cigar Overlaps::get_overlap(const HandleGraph& graph, handle_t a, handle_t b) const {
    Cigar cigar;
    auto it = overlaps.find(graph.edge_handle(a, b));
    if (it != overlaps.end()) {
        if (it->first.first == a) {
            cigar = it->second;
        }
        else {
            cigar = it->second.reverse();
        }
    }
    return cigar;
}


