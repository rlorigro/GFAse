#ifndef GFASE_HAMILTONIAN_CHAINER_HPP
#define GFASE_HAMILTONIAN_CHAINER_HPP

#include "MultiContactGraph.hpp"
#include "IncrementalIdMap.hpp"
#include "hash_graph.hpp"

using bdsg::MutablePathHandleGraph;
using bdsg::HandleGraph;
using handlegraph::nid_t;
using handlegraph::path_handle_t;

#include <unordered_set>
#include <utility>
#include <array>
#include <vector>

using std::unordered_set;
using std::pair;
using std::vector;
using std::array;

namespace gfase {

class HamiltonianChainer {
public:
    
    HamiltonianChainer() = default;
    ~HamiltonianChainer() = default;
    
    // expects nodes that do not have alts to have been deleted from
    // the contact graph
    void generate_chain_paths(MutablePathHandleGraph& graph,
                              const MultiContactGraph& contact_graph) const;
    
    // the number of iterations we allow in a Hamiltonian path problem
    size_t hamiltonian_max_iters = 5000;
    
private:
    
    // returns the confident left and right sides of the allelic walk through a component.
    // if there is a full-length walk, only the first vector in the pair is filled.
    // walks are all oriented away from the starts, toward the ends.
    // also records whether the alleles were generated with a hamiltonian walk
    pair<vector<handle_t>, vector<handle_t>> generate_allelic_semiwalks(const HandleGraph& graph,
                                                                        const unordered_set<nid_t>& in_phase_nodes,
                                                                        const unordered_set<nid_t>& out_phase_nodes,
                                                                        const unordered_set<handle_t>& starts,
                                                                        const unordered_set<handle_t>& ends,
                                                                        bool& resolved_hamiltonian) const;
    
    // generate the path names we use to mark haplotypes
    static string phase_path_name(int haplotype, int path_id);
    
    // the maximum length walk using only the allowed nodes before encountering
    // a branch to multiple allowed nodes
    vector<handle_t> max_unambiguous_path(const HandleGraph& graph, handle_t start,
                                          const unordered_set<nid_t>& allowed_nodes) const;
    
    // find the longest non-branching path that consists of only unphased (diploid) nodes starting at the handle
    vector<handle_t> walk_diploid_unipath(const HandleGraph& graph, const MultiContactGraph& contact_graph,
                                          handle_t handle, bool go_left) const;
    
    // split phase paths into 2 at self-loops
    void break_self_looping_phase_paths(MutablePathHandleGraph& graph,
                                        array<unordered_set<path_handle_t>, 2>& phase_paths,
                                        array<int, 2>& next_path_ids) const;
    
    // check if we can extend phase path unambiguously into paths into components that
    // we couldn't previously phase
    void extend_unambiguous_phase_paths(MutablePathHandleGraph& graph,
                                        array<unordered_set<path_handle_t>, 2>& phase_paths,
                                        array<int, 2>& next_path_ids) const;
};

}


#endif //GFASE_HAMILTONIAN_CHAINER_HPP
