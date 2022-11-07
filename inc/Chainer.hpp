#ifndef GFASE_CHAINER_HPP
#define GFASE_CHAINER_HPP

#include "MultiContactGraph.hpp"
#include "IncrementalIdMap.hpp"
#include "hash_graph.hpp"
#include "Filesystem.hpp"

#include "hash_graph.hpp"

using bdsg::MutablePathDeletableHandleGraph;
using bdsg::HashGraph;
using bdsg::handle_t;
using bdsg::nid_t;

using ghc::filesystem::path;

#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <stdexcept>
#include <fstream>
#include <utility>
#include <limits>
#include <vector>
#include <array>
#include <set>
#include <deque>

using std::numeric_limits;
using std::unordered_set;
using std::unordered_map;
using std::runtime_error;
using std::exception;
using std::to_string;
using std::function;
using std::ifstream;
using std::ofstream;
using std::string;
using std::vector;
using std::array;
using std::pair;
using std::cerr;
using std::ref;
using std::set;
using std::deque;

namespace gfase {


using chain_t = deque <set <nid_t> >;


class Chainer {
    unordered_map<nid_t,nid_t> node_pairs;
    unordered_set<nid_t> diploid_nodes;
    unordered_set<nid_t> diploid_tip_nodes;
    unordered_set<nid_t> haploid_nodes;

public:
    // Constructor
    Chainer()=default;

    // Chaining
    void find_chainable_nodes(const HandleGraph& graph, const IncrementalIdMap<string>& id_map);

    void get_chain(const HandleGraph& graph,
                   const nid_t& start_node,
                   deque <set <nid_t> >& chain);

    void validate_chain_paths();

    void for_each_chain(
            HandleGraph& graph,
            const function<void(chain_t& chain)>& f);

    void generate_chain_paths(
            MutablePathDeletableHandleGraph& graph,
            const IncrementalIdMap<string>& id_map,
            const MultiContactGraph& contact_graph);

    // IO
    void write_chainable_nodes_to_bandage_csv(path output_dir, const IncrementalIdMap<string>& id_map) const;
};

}


#endif //GFASE_CHAINER_HPP
