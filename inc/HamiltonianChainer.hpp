#ifndef GFASE_HAMILTONIAN_CHAINER_HPP
#define GFASE_HAMILTONIAN_CHAINER_HPP

#include "MultiContactGraph.hpp"
#include "IncrementalIdMap.hpp"
#include "hash_graph.hpp"

using bdsg::MutablePathDeletableHandleGraph;
using bdsg::HandleGraph;
using handlegraph::nid_t;

#include <unordered_set>
#include <utility>

using std::unordered_set;
using std::pair;
using std::vector;

namespace gfase {

class HamiltonianChainer {
public:
    
    HamiltonianChainer() = default;
    ~HamiltonianChainer() = default;
    
    // expects nodes that do not have alts to have been deleted from
    // the contact graph
    void generate_chain_paths(MutablePathDeletableHandleGraph& graph,
                              const IncrementalIdMap<string>& id_map,
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
                                                                        const unordered_set<handle_t>& ends
                                                                        bool& resolved_hamiltonian) const;
    
    // the maximum length walk using only the allowed nodes before encountering
    // a branch to multiple allowed nodes
    vector<handle_t> max_unambiguous_path(const HandleGraph& graph, handle_t start,
                                          const unordered_set<nid_t>& allowed_nodes) const;
    
    vector<handle_t> walk_diploid_unipath(const HandleGraph& graph, const MultiContactGraph& contact_graph
                                          handle_t handle, bool go_left) const;
    
};

}


#endif //GFASE_HAMILTONIAN_CHAINER_HPP
