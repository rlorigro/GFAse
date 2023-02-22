#ifndef GFASE_HAMILTONIAN_CHAINER_HPP
#define GFASE_HAMILTONIAN_CHAINER_HPP

#include "MultiContactGraph.hpp"
#include "IncrementalIdMap.hpp"
#include "Chainer.hpp"
#include "hash_graph.hpp"

using bdsg::MutablePathDeletableHandleGraph;
using bdsg::HandleGraph;
using handlegraph::nid_t;
using handlegraph::path_handle_t;

#include <unordered_set>
#include <utility>
#include <array>
#include <vector>
#include <string>

using std::unordered_set;
using std::pair;
using std::vector;
using std::array;
using std::string;

namespace gfase {

class HamiltonianChainer : public AbstractChainer {
public:
    
    HamiltonianChainer() = default;
    ~HamiltonianChainer() = default;
    
    // add phased haplotype paths to the graph
    void generate_chain_paths(MutablePathDeletableHandleGraph& graph,
                              const IncrementalIdMap<string>& id_map, // this isn't used, but it keeps a consistent interface
                              const MultiContactGraph& contact_graph);
    
    // is this the name of one of the phased haplotype paths?
    bool has_phase_chain(const string& name) const;
    // what is the partition of the phased hapotype path?
    int8_t get_partition(const string& name) const;
    // after generating chain paths, write them to bandage annotations
    void write_chaining_results_to_bandage_csv(path output_dir, const IncrementalIdMap<string>& id_map,
                                               const MultiContactGraph& contact_graph) const;
    
    // the number of iterations we allow in a Hamiltonian path problem
    size_t hamiltonian_max_iters = 5000;
    
    // we will discard putative haplotypes consisting of less than this proportion of
    // phased haploid sequence
    double min_haploid_proportion = 0.1;
    
private:
    
    // keep track of which of the paths are the phase paths we added
    array<unordered_set<path_handle_t>, 2> phase_paths;
    // we remember the graph so that we can use some of its in-built indexes in the public inferface
    // TODO: ugly solution
    const PathHandleGraph* path_graph = nullptr;
    
    // returns the confident left and right sides of the allelic walk through a component.
    // if there is a full-length walk, only the first vector in the pair is filled.
    // walks are all oriented away from the starts, toward the ends.
    // also records whether the alleles were generated with a hamiltonian walk
    pair<vector<handle_t>, vector<handle_t>> generate_allelic_semiwalks(const HandleGraph& graph,
                                                                        const IncrementalIdMap<string>& id_map,
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
    
    // check if we can extend phase path unambiguously into paths into components that
    // we couldn't previously phase
    void extend_unambiguous_phase_paths(MutablePathDeletableHandleGraph& graph,
                                        const MultiContactGraph& contact_graph);
    
    // split phase paths into 2 at self-loops
    void break_self_looping_phase_paths(MutablePathDeletableHandleGraph& graph,
                                        array<int, 2> next_path_ids);
    
    // remove phase paths that don't actually have any phased, haploid content
    void purge_null_phase_paths(MutablePathDeletableHandleGraph& graph,
                                const MultiContactGraph& contact_graph);
    
    // remove phase paths that don't meet the minimum phased proportion, unless they are extensions of
    // alleles that exist but weren't fully resolved, or are the alternate allele on a bridge for
    // a path that does meet the minimum
    void purge_mostly_unphased_phase_paths(MutablePathDeletableHandleGraph& graph,
                                           const MultiContactGraph& contact_graph,
                                           const unordered_map<path_handle_t, unordered_set<path_handle_t>>>& bridge_overlap_graph,
                                           const unordered_map<path_handle_t, unordered_set<path_handle_t>>>& broken_allele_graph);
};

}


#endif //GFASE_HAMILTONIAN_CHAINER_HPP
