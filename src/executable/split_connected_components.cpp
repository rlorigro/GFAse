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
using gfase::write_connected_components_to_gfas;
using gfase::handle_graph_to_gfa;
using gfase::print_graph_paths;
using gfase::plot_graph;
using ghc::filesystem::path;
using ghc::filesystem::create_directories;

using bdsg::HashGraph;
using bdsg::MutablePathMutableHandleGraph;
using bdsg::MutablePathDeletableHandleGraph;
using handlegraph::path_handle_t;
using handlegraph::step_handle_t;
using handlegraph::handle_t;

using std::string;
using std::cout;
using std::cerr;


// TODO: fix for whole genome, missing nodes/links in path ??


void split_gfa_components(path gfa_path){
    HashGraph graph;
    IncrementalIdMap<string> id_map;
    Overlaps overlaps(graph);

    cerr << "Loading GFA..." << '\n';
    gfa_to_handle_graph(graph, id_map, overlaps, gfa_path);

    path output_directory = gfa_path.parent_path() / gfa_path.stem();

    cerr << "Writing subgraph GFAs to: " << output_directory << '\n';
    create_directories(output_directory);

    write_connected_components_to_gfas(graph, id_map, output_directory);
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

    split_gfa_components(gfa_path);

    return 0;
}
