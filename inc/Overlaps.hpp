#ifndef GFASE_OVERLAPS_HPP
#define GFASE_OVERLAPS_HPP

#include "bdsg/hash_graph.hpp"
#include <unordered_map>
#include <string>
#include <vector>
#include <cstdint>
#include <array>


using handlegraph::handle_t;
using handlegraph::edge_t;
using std::unordered_map;
using bdsg::HandleGraph;
using std::string;
using std::vector;
using std::array;
using std::pair;


class Cigar; // forward declaration

/*
 * Single operation in a CIGAR string
 */
class CigarOperation{
public:
    CigarOperation(uint32_t length, char type);
    CigarOperation() = default;
    ~CigarOperation() = default;
    
    char type() const;
    uint32_t length() const;
    
private:
    
    // Map from all possible cigar chars to 0-8
    static const array<uint8_t,128> cigar_code;
    
    // Map back to the human readable characters
    static const array<char,9> cigar_type;
    
    // Check if a cigar operation consumes a base in the reference or query
    static const array<bool,9> is_ref_move;
    static const array<bool,9> is_query_move;
    
    /// Attributes ///
    uint8_t code;
    uint32_t op_length;
};

/*
 * A parsed CIGAR string
 */
class Cigar {
public:
    // parse a CIGAR string
    Cigar(const std::string& cigar_string);
    Cigar() = default;
    ~Cigar() = default;
    
    // get the CIGAR string for this alignment
    std::string get_string() const;
    
    // return the CIGAR string for the same alignment, but reversed
    // and with query/ref roles swapped
    Cigar reverse() const;
    
    // number of operations
    size_t size() const;
    // i-th operation
    const CigarOperation& at(size_t i) const;
    
    // foreach interface
    using iterator = std::vector<CigarOperation>::const_iterator;
    iterator begin() const;
    iterator end() const;
    
private:
    
    std::vector<CigarOperation> operations;
    
};

/*
 * Records overlaps between nodes of a graph
 */
class Overlaps {
public:
    Overlaps(const HandleGraph& graph);
    Overlaps() = default;
    ~Overlaps() = default;
    
    
    void record_overlap(handle_t a, handle_t b, const string& cigar);
    
    // true if there are no overlaps recorded
    bool is_blunt() const;
    
    // true if there is an edge from a to b with an overlap
    bool has_overlap(handle_t a, handle_t b) const;
    
    // get the (oriented) overlap of a onto b
    Cigar get_overlap(handle_t a, handle_t b) const;
    
private:
    
    const HandleGraph* graph;
    unordered_map<edge_t, Cigar> overlaps;
};


#endif //GFASE_OVERLAPS_HPP
