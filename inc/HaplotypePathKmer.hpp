#ifndef GFASE_HAPLOTYPEPATHKMER_HPP
#define GFASE_HAPLOTYPEPATHKMER_HPP

#include "bdsg/hash_graph.hpp"

#include <deque>

using bdsg::HashGraph;
using handlegraph::MutablePathMutableHandleGraph;
using handlegraph::MutablePathDeletableHandleGraph;
using handlegraph::PathHandleGraph;
using handlegraph::path_handle_t;
using handlegraph::step_handle_t;
using handlegraph::handle_t;

using std::deque;


namespace gfase {


bool is_haplotype_bubble(const PathHandleGraph& graph, step_handle_t s);


class HaplotypePathKmer {
public:
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

    /// Methods ///
    HaplotypePathKmer(const PathHandleGraph& graph, const path_handle_t& path, size_t k);

    // skip - jump to a position (step_handle_t, size_t) in the path and update internal records
    void initialize(step_handle_t s, size_t index);

    // step - walk a single bp forward and update internal records
    bool step();

    bool update_has_diploid();
};


}

#endif //GFASE_HAPLOTYPEPATHKMER_HPP
