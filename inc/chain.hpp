#ifndef GFASE_CHAIN_HPP
#define GFASE_CHAIN_HPP

#include "IncrementalIdMap.hpp"
#include "gfa_to_handle.hpp"
#include "graph_utility.hpp"
#include "BubbleGraph.hpp"
#include "Bipartition.hpp"
#include "hash_graph.hpp"
#include "Filesystem.hpp"
#include "Timer.hpp"
#include "CLI11.hpp"
#include "Sam.hpp"
#include "Bam.hpp"

using gfase::for_element_in_sam_file;
using gfase::unpaired_mappings_t;
using gfase::paired_mappings_t;
using gfase::contact_map_t;

using gfase::gfa_to_handle_graph;
using gfase::IncrementalIdMap;
using gfase::Bipartition;
using gfase::BubbleGraph;
using gfase::SamElement;
using gfase::Bubble;
using gfase::Timer;
using gfase::Bam;

using bdsg::HashGraph;


using ghc::filesystem::path;
using CLI::App;

#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <stdexcept>
#include <fstream>
#include <utility>
#include <atomic>
#include <thread>
#include <limits>
#include <bitset>
#include <vector>
#include <mutex>
#include <array>
#include <set>

using std::numeric_limits;
using std::unordered_set;
using std::unordered_map;
using std::runtime_error;
using std::streamsize;
using std::exception;
using std::to_string;
using std::function;
using std::ifstream;
using std::ofstream;
using std::shuffle;
using std::string;
using std::vector;
using std::bitset;
using std::thread;
using std::atomic;
using std::array;
using std::mutex;
using std::pair;
using std::stoi;
using std::cerr;
using std::cref;
using std::ref;
using std::set;


namespace gfase{

void generate_ploidy_criteria_from_bubble_graph(
        const PathHandleGraph& graph,
        const IncrementalIdMap<string>& id_map,
        const BubbleGraph& bubble_graph,
        unordered_set<nid_t>& diploid_nodes);


void merge_diploid_singletons(const BubbleGraph& bubble_graph, Bipartition& chain_bipartition);


void write_chaining_info_to_file(
        path output_dir,
        const Bipartition& ploidy_bipartition,
        const Bipartition& chain_bipartition,
        const PathHandleGraph& graph,
        const IncrementalIdMap<string>& id_map,
        const string& filename_prefix,
        size_t component_index,
        bool write_gfa = true);


void chain_phased_gfa(
        MutablePathDeletableHandleGraph& graph,
        IncrementalIdMap<string>& id_map,
        Overlaps& overlaps,
        const BubbleGraph& bubble_graph,
        path output_dir,
        bool write_gfa = true,
        bool write_fasta  = true);

}

#endif //GFASE_CHAIN_HPP
