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


void unzip_gfa(path gfa_path){
    HashGraph graph;
    IncrementalIdMap<string> id_map;
    Overlaps overlaps;

    gfa_to_handle_graph(graph, id_map, overlaps, gfa_path);

    // Output an image of the graph, can be uncommented for debugging
//    plot_graph(graph, "test_unzip_unedited");

    vector<HashGraph> connected_components;
    vector <IncrementalIdMap<string> > connected_component_ids;
    vector<Overlaps> connected_component_overlaps;

    split_connected_components(graph, id_map, overlaps, connected_components, connected_component_ids, connected_component_overlaps);

    for (size_t i=0; i<connected_components.size(); i++){
        cerr << "Component " << to_string(i) << '\n';
        print_graph_paths(connected_components[i], connected_component_ids[i]);

        unzip(connected_components[i], connected_component_ids[i], connected_component_overlaps[i]);

        string filename_prefix = "component_" + to_string(i) + "_unzipped";
        ofstream file(filename_prefix + ".gfa");
        handle_graph_to_gfa(connected_components[i], connected_component_ids[i], file);
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

    unzip_gfa(gfa_path);

    return 0;
}
