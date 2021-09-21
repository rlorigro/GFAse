#ifndef GFASE_HAPLOTYPEPATHKMER_HPP
#define GFASE_HAPLOTYPEPATHKMER_HPP

#include "bdsg/hash_graph.hpp"

#include <functional>
#include <deque>
#include <string>

using bdsg::HashGraph;
using handlegraph::MutablePathMutableHandleGraph;
using handlegraph::MutablePathDeletableHandleGraph;
using handlegraph::PathHandleGraph;
using handlegraph::path_handle_t;
using handlegraph::step_handle_t;
using handlegraph::handle_t;

using std::function;
using std::deque;
using std::string;


namespace gfase {


bool is_haplotype_bubble(const PathHandleGraph& graph, step_handle_t s);


class HaplotypePathKmer {
private:
    /// Attributes ///
    const PathHandleGraph& graph;
    deque <step_handle_t> steps;
    deque <size_t> lengths;
    deque<bool> is_diploid;
    bool has_diploid;
    path_handle_t path;

    size_t start_index;
    size_t stop_index;

    step_handle_t terminal_step;
    size_t terminal_index;

    // K can be any size, so bit shift operations are not a simple option anymore
    deque<char> sequence;
    size_t k;

public:
    /// Methods ///
    HaplotypePathKmer(const PathHandleGraph& graph, const path_handle_t& path, size_t k);

    // skip - jump to a position (step_handle_t, size_t) in the path and update internal records
    void initialize(step_handle_t s, size_t index);

    // step - walk a single bp forward and update internal records
    bool step();

    void for_each_haploid_kmer(const function<void(deque<char>& sequence)>& f);

    bool update_has_diploid();

    step_handle_t get_step_of_kmer_start() const;
    size_t get_index_of_kmer_start() const;
    void print();
};


}

#endif //GFASE_HAPLOTYPEPATHKMER_HPP
