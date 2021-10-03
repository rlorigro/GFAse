#include "HaplotypePathKmer.hpp"
#include "IncrementalIdMap.hpp"
#include "gfa_to_handle.hpp"
#include "handle_to_gfa.hpp"
#include "graph_utility.hpp"
#include "Filesystem.hpp"
#include "GfaReader.hpp"
#include "CLI11.hpp"

#include "bdsg/hash_graph.hpp"

#include <string>

using gfase::HaplotypePathKmer;
using gfase::IncrementalIdMap;
using gfase::find_singleton_adjacent_handle;
using gfase::for_each_connected_component;
using gfase::split_connected_components;
using gfase::handle_graph_to_gfa;
using gfase::find_diploid_paths;
using gfase::print_graph_paths;
using gfase::extend_paths;
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


void extract_haplotype_kmers_from_gfa(path gfa_path, size_t k){
    HashGraph graph;
    IncrementalIdMap<string> id_map;

    gfa_to_handle_graph(graph, id_map, gfa_path);

    plot_graph(graph, "start_graph");

    cerr << "Identifying diploid paths..." << '\n';

    unordered_set<string> diploid_path_names;
    find_diploid_paths(graph, diploid_path_names);

    cerr << "Extending paths into haploid regions..." << '\n';

    vector<pair<path_handle_t, handle_t> > to_be_prepended;
    vector<pair<path_handle_t, handle_t> > to_be_appended;

    extend_paths(graph, to_be_prepended, to_be_appended);

    cerr << "Iterating path kmers..." << '\n';

    // Iterate paths and for each node, collect kmers if node is only covered by one path
    for (auto& path_name: diploid_path_names) {
        auto p = graph.get_path_handle(path_name);

        cerr << ">" << path_name << '\n';

        HaplotypePathKmer kmer(graph, p, k);

        kmer.for_each_haploid_kmer([&](const deque<char>& sequence){
            for (const auto& c: sequence){
                cerr << c;
            }
            cerr << '\n';
        });
    }
}


int main (int argc, char* argv[]){
    path gfa_path;
    size_t k;

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

    CLI11_PARSE(app, argc, argv);

    extract_haplotype_kmers_from_gfa(gfa_path, k);

    return 0;
}
