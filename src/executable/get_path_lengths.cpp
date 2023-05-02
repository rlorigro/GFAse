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
using handlegraph::path_handle_t;
using handlegraph::step_handle_t;
using handlegraph::handle_t;

using std::string;
using std::cout;
using std::cerr;


int get_path_lengths(path gfa_path){
    HashGraph graph;
    IncrementalIdMap<string> id_map;
    Overlaps overlaps;

    gfa_to_handle_graph(graph, id_map, overlaps, gfa_path, false);

    graph.for_each_path_handle([&](const path_handle_t& p) {
        cerr << graph.get_path_name(p) << ": ";
        size_t l = 0;
        graph.for_each_step_in_path(p, [&](const step_handle_t& s){
            auto h = graph.get_handle_of_step(s);

            l += graph.get_length(h);
        });

        cerr << l << '\n';
    });


    return 0;
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

    get_path_lengths(gfa_path);

    return 0;
}
