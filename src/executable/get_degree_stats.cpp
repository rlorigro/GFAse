#include "IncrementalIdMap.hpp"
#include "gfa_to_handle.hpp"
#include "handle_to_gfa.hpp"
#include "graph_utility.hpp"
#include "Filesystem.hpp"
#include "GfaReader.hpp"
#include "CLI11.hpp"

#include "bdsg/hash_graph.hpp"

#include <string>

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
using bdsg::path_handle_t;
using bdsg::step_handle_t;
using bdsg::handle_t;

using std::string;
using std::cout;
using std::cerr;

void get_degree_stats(path gfa_path){
    HashGraph graph;
    IncrementalIdMap<string> id_map;

    gfa_to_handle_graph(graph, id_map, gfa_path, false);

    map<size_t,size_t> degree_counts;

    graph.for_each_handle([&](const handle_t& h){
        size_t d;
        d = graph.get_degree(h, true);
        degree_counts[d]++;
        d = graph.get_degree(h, false);
        degree_counts[d]++;
    });

    for (auto& [d,n]: degree_counts){
        cout << d << '\t' << n << '\n';
    }
}


int main (int argc, char* argv[]){
    path gfa_path;

    CLI::App app{"App description"};

    app.add_option(
                    "-i,--input_gfa",
                    gfa_path,
                    "Path to GFA containing phased non-overlapping segments")
            ->required();

    CLI11_PARSE(app, argc, argv);

    get_degree_stats(gfa_path);

    return 0;
}
