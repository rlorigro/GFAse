#include "IncrementalIdMap.hpp"
#include "gfa_to_handle.hpp"
#include "handle_to_gfa.hpp"
#include "graph_utility.hpp"
#include "BubbleGraph.hpp"
#include "Filesystem.hpp"
#include "CLI11.hpp"

#include "bdsg/hash_graph.hpp"

#include <string>

using gfase::IncrementalIdMap;
using gfase::BubbleGraph;
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


void find_bubbles_in_gfa(path output_dir, path gfa_path){
    if (exists(output_dir)){
        throw runtime_error("ERROR: output directory exists already");
    }
    else {
        create_directories(output_dir);
    }

    HashGraph graph;
    IncrementalIdMap<string> id_map;
    Overlaps overlaps;

    gfa_to_handle_graph(graph, id_map, overlaps, gfa_path, false);

    BubbleGraph bubble_graph(graph);

    path output_path = output_dir / "bubbles.csv";
    bubble_graph.write_bandage_csv(output_path, id_map);
}


int main (int argc, char* argv[]){
    path gfa_path;
    path output_dir;

    CLI::App app{"App description"};

    app.add_option(
            "-g,--gfa",
            gfa_path,
            "Path to GFA containing assembly graph to be phased")
            ->required();

    app.add_option(
            "-o,--output_dir",
            output_dir,
            "Path to (nonexistent) directory where output will be stored")
            ->required();


    CLI11_PARSE(app, argc, argv);

    find_bubbles_in_gfa(output_dir, gfa_path);

    return 0;
}
