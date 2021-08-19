#include "IncrementalIdMap.hpp"
#include "gfa_to_handle.hpp"
#include "handle_to_gfa.hpp"
#include "Filesystem.hpp"
#include "GfaReader.hpp"
#include "CLI11.hpp"

#include "bdsg/hash_graph.hpp"

#include <string>

using gfase::IncrementalIdMap;
using gfase::handle_graph_to_gfa;
using ghc::filesystem::path;

using bdsg::HashGraph;
using bdsg::MutablePathMutableHandleGraph;
using handlegraph::path_handle_t;
using handlegraph::step_handle_t;
using handlegraph::handle_t;

using std::string;
using std::cout;
using std::cerr;


void run_command(string& argument_string){
    int exit_code = system(argument_string.c_str());

    if (exit_code != 0){
        throw runtime_error("ERROR: command failed to run: " + argument_string);
    }
}


void plot_graph(HandleGraph& graph, string filename_prefix){
    ofstream test_output(filename_prefix + ".gfa");
    handle_graph_to_gfa(graph, test_output);
    test_output.close();

    if (graph.get_node_count() < 200) {
        string command = "vg convert -g " + filename_prefix + ".gfa -p | vg view -d - | dot -Tpng -o "
                         + filename_prefix + ".png";

        cerr << "Running: " << command << '\n';

        run_command(command);
    }
}


void split_into_connected_components(HashGraph& graph, IncrementalIdMap<string> id_map){


}


void unzip(path gfa_path){
    HashGraph graph;
    IncrementalIdMap<string> id_map;

    gfa_to_handle_graph(graph, id_map, gfa_path);

    // Output an image of the graph, can be uncommented for debugging
    plot_graph(graph, "test_unzip_unedited");

    unordered_set<handle_t> to_be_destroyed;

    cout << graph.get_path_count() << '\n';
    vector<string> path_names;

    graph.for_each_path_handle([&](const path_handle_t& p) {
        string path_sequence;

        graph.for_each_step_in_path(p, [&](const step_handle_t s){
            handle_t h = graph.get_handle_of_step(s);
            nid_t n = graph.get_id(h);

            string sequence = graph.get_sequence(h);
            path_sequence += sequence;

            string name = id_map.get_name(n);
            cout << '\t' << name << " " << sequence << '\n';

            to_be_destroyed.emplace(h);
        });

        // TODO: track provenance?
        handle_t haplotype_handle = graph.create_handle(path_sequence);
        auto path_start_handle = graph.get_handle_of_step(graph.path_begin(p));
        auto path_stop_handle = graph.get_handle_of_step(graph.path_back(p));

        // Find neighboring nodes for the path and create edges to the new haplotype node (LEFT)
        graph.follow_edges(path_start_handle, true, [&](const handle_t& other){
            graph.create_edge(other, haplotype_handle);
        });

        // Find neighboring nodes for the path and create edges to the new haplotype node (RIGHT)
        graph.follow_edges(path_stop_handle, false, [&](const handle_t& other){
            graph.create_edge(haplotype_handle, other);
        });
    });

    // Output an image of the graph, can be uncommented for debugging
    plot_graph(graph, "test_unzip_duplicated");

    // Destroy all the paths involved in the unzipping
    graph.for_each_path_handle([&](const path_handle_t& p) {
        graph.destroy_path(p);
    });

    // Destroy the nodes that were duplicated into haplotypes
    for (auto& handle: to_be_destroyed){
        graph.destroy_handle(handle);
    }

    // Output an image of the graph, can be uncommented for debugging
    plot_graph(graph, "test_unzip_final");
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

    unzip(gfa_path);

    return 0;
}
