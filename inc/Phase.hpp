#ifndef GFASE_PHASE_HPP
#define GFASE_PHASE_HPP

#include "HaplotypePathKmer.hpp"
#include "IncrementalIdMap.hpp"
#include "gfa_to_handle.hpp"
#include "handle_to_gfa.hpp"
#include "graph_utility.hpp"
#include "Filesystem.hpp"
#include "GfaReader.hpp"
#include "KmerSets.hpp"
#include "CLI11.hpp"

#include "bdsg/hash_graph.hpp"

#include <string>

using ghc::filesystem::path;

using bdsg::HashGraph;
using bdsg::MutablePathMutableHandleGraph;
using bdsg::MutablePathDeletableHandleGraph;
using handlegraph::path_handle_t;
using handlegraph::step_handle_t;
using handlegraph::handle_t;

using std::string;
using std::cout;
using std::cerr;


namespace gfase {

// Find any nodes that are adjacent to the beginning and end of a path, as long as they are the only adjacent node
pair<handle_t, bool> find_singleton_adjacent_handle(const PathHandleGraph& graph, const handle_t& h, bool left);


void extend_paths(MutablePathMutableHandleGraph& graph);


void phase_haplotype_paths(path gfa_path, size_t k, path paternal_kmers, path maternal_kmers, char path_delimiter='.');

}

#endif //GFASE_PHASE_HPP
