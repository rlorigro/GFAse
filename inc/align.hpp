#ifndef GFASE_ALIGN_HPP
#define GFASE_ALIGN_HPP

#include "IncrementalIdMap.hpp"
#include "gfa_to_handle.hpp"
#include "graph_utility.hpp"
#include "ContactGraph.hpp"
#include "Bipartition.hpp"
#include "hash_graph.hpp"
#include "Filesystem.hpp"
#include "Sequence.hpp"
#include "Hasher2.hpp"
#include "Timer.hpp"
#include "CLI11.hpp"
#include "Sam.hpp"
#include "Bam.hpp"
#include "minimap.h"
#include "SvgPlot.hpp"

using ghc::filesystem::path;
using CLI::App;

#include <unordered_map>

using std::unordered_map;

namespace gfase {


void print_minimap_alignment_block(
        mm_mapopt_t& map_options,
        mm_idx_t* mi,
        mm_reg1_t* r2,
        const string& name,
        const string& query);


void map_sequence_pair(
        const string& target_name,
        const string& target_sequence,
        const string& query_name,
        const string& query_sequence,
        AlignmentChain& result);


void construct_alignment_graph(
        const vector<HashResult>& to_be_aligned,
        const vector<Sequence>& sequences,
        const IncrementalIdMap<string>& id_map,
        ContactGraph& alignment_graph,
        double min_similarity,
        mutex& output_mutex,
        atomic<size_t>& global_index);


void get_alignment_candidates(
        const vector<Sequence>& sequences,
        const IncrementalIdMap<string>& id_map,
        vector <HashResult>& to_be_aligned,
        path output_dir,
        size_t n_threads,
        double sample_rate,
        size_t k,
        size_t n_iterations,
        size_t max_hits,
        double min_ab_over_a,
        double min_ab_over_b
        );


void get_best_overlaps(
        double min_similarity,
        const IncrementalIdMap<string>& id_map,
        ContactGraph& alignment_graph,
        ContactGraph& symmetrical_alignment_graph
        );


void write_alignment_results_to_file(
        const IncrementalIdMap<string>& id_map,
        const ContactGraph& alignment_graph,
        const ContactGraph& symmetrical_alignment_graph,
        path output_dir);


}

#endif //GFASE_ALIGN_HPP
