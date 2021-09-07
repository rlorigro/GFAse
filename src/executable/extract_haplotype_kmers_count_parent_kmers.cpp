#include "HaplotypePathKmer.hpp"
#include "IncrementalIdMap.hpp"
#include "gfa_to_handle.hpp"
#include "handle_to_gfa.hpp"
#include "GraphUtility.hpp"
#include "Filesystem.hpp"
#include "GfaReader.hpp"
#include "kmer_unordered_set.hpp"
#include "CLI11.hpp"

#include "bdsg/hash_graph.hpp"

#include <string>

using gfase::HaplotypePathKmer;
using gfase::IncrementalIdMap;
using gfase::for_each_connected_component;
using gfase::split_connected_components;
using gfase::handle_graph_to_gfa;
using gfase::print_graph_paths;
using gfase::plot_graph;
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


// Find any nodes that are adjacent to the beginning and end of a path, as long as they are the only adjacent node
pair<handle_t, bool> find_singleton_adjacent_handle(const PathHandleGraph& graph, const handle_t& h, bool left){
    handle_t adjacent_handle;
    size_t n_adjacent = 0;
    bool success = false;

    graph.follow_edges(h, left, [&](const handle_t& other_handle){
        if (n_adjacent > 0){
            return false;
        }

        adjacent_handle = other_handle;
        n_adjacent++;

        return true;
    });

    // Only return true if there was exactly one adjacent handle
    if (n_adjacent == 1){
        success = true;
    }

    return {adjacent_handle, success};
}


void extend_paths(MutablePathMutableHandleGraph& graph){
    vector<pair <path_handle_t, handle_t> > to_be_prepended;
    vector<pair <path_handle_t, handle_t> > to_be_appended;

    graph.for_each_path_handle([&](const path_handle_t& p){
        auto begin_handle = graph.get_handle_of_step(graph.path_begin(p));
        auto end_handle = graph.get_handle_of_step(graph.path_back(p));

        auto left_result = find_singleton_adjacent_handle(graph, begin_handle, true);
        auto right_result = find_singleton_adjacent_handle(graph, end_handle, false);

        if (left_result.second){
            size_t n_paths = 0;
            graph.for_each_step_on_handle(left_result.first, [&](const step_handle_t& s){
                n_paths++;
            });

            // Verify that there is not already some other path covering the adjacent node
            if (n_paths == 0) {
                to_be_prepended.emplace_back(p, left_result.first);
            }
            else{
                throw runtime_error("ERROR: node left of haplotype path " + graph.get_path_name(p) + " is covered by other paths");
            }
        }

        if (right_result.second){
            size_t n_paths = 0;
            graph.for_each_step_on_handle(right_result.first, [&](const step_handle_t& s){
                n_paths++;
            });

            // Verify that there is not already some other path covering the adjacent node
            if (n_paths == 0) {
                to_be_appended.emplace_back(p, right_result.first);
            }
            else{
                throw runtime_error("ERROR: node right haplotype path " + graph.get_path_name(p) + " is covered by other paths");
            }
        }
    });

    for (auto& item: to_be_prepended){
        // PREpend the LEFT side node if it meets the conditions
        graph.prepend_step(item.first, item.second);
    }

    for (auto& item: to_be_appended){
        // Append the RIGHT side node if it meets the conditions
        graph.append_step(item.first, item.second);
    }
}


void extract_haplotype_kmers_from_gfa(path gfa_path, size_t k, path paternal_kmers, path maternal_kmers){
    HashGraph graph;
    IncrementalIdMap<string> id_map;
    KmerSets ks(paternal_kmers, maternal_kmers);

    gfa_to_handle_graph(graph, id_map, gfa_path);

    // plot_graph(graph, "start_graph");

    cerr << "Extending paths into haploid regions..." << '\n';

    extend_paths(graph);

    cerr << "Iterating path kmers..." << '\n';

    bool prev_has_diploid;

    // Iterate paths and for each node, collect kmers if node is only covered by one path
    graph.for_each_path_handle([&](const path_handle_t& p){
        cerr << ">" << graph.get_path_name(p) << '\n';
 
        HaplotypePathKmer kmer(graph, p, k);

        while (kmer.step()){
            if (kmer.has_diploid){
                string kmer_string;
                for (auto& c: kmer.sequence){
                    cerr << c; 
                    kmer_string+=c;
                }
                // compare kmer to parental kmers
                ks.find_haplotype_single_kmer_count(kmer_string);
                cerr << '\n'; 

            }
            else{
                auto h = graph.get_handle_of_step(kmer.steps.back());
                auto length = graph.get_length(h);

                if (length > k + 1){
                    kmer.initialize(kmer.steps.back(), length - k + 1);
                    prev_has_diploid = false;
                }
            }

            prev_has_diploid = kmer.has_diploid;
        }
    });
}


int main (int argc, char* argv[]){
    path gfa_path;
    size_t k;
    path paternal_kmers;
    path maternal_kmers;

    CLI::App app{"App description"};

    app.add_option(
            "-i,--input_gfa",
            gfa_path,
            "Path to GFA containing phased non-overlapping segments")
            ->required();

    app.add_option(
            "-k,--kmer_size",
            k,
            "Length of kmer (k) to use")
            ->required();

    app.add_option(
            "-p,--paternal_kmer_fa",
            paternal_kmers,
            "paternal kmers in .fa format")
            ->required();

    app.add_option(
            "-m,--maternal_kmer_fa",
            maternal_kmers,
            "maternal kmers in .fa format")
            ->required();

    CLI11_PARSE(app, argc, argv);

    extract_haplotype_kmers_from_gfa(gfa_path, k, paternal_kmers, maternal_kmers);

    return 0;
}
