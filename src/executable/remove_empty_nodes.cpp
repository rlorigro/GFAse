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


void splice_node_neighbors(MutableHandleGraph& graph, handle_t h){
    vector <handle_t> in_nodes;
    vector <handle_t> out_nodes;

    graph.follow_edges(h, true, [&](const handle_t& h_other){
        in_nodes.emplace_back(h_other);
    });

    graph.follow_edges(h, false, [&](const handle_t& h_other){
        out_nodes.emplace_back(h_other);
    });

    for (auto& h_in: in_nodes){
        for (auto& h_out: out_nodes){
            graph.create_edge(h_in, h_out);
        }
    }
}


void splice_out_empty_nodes(MutablePathDeletableHandleGraph& graph){
    vector <path_handle_t> paths;

    graph.for_each_path_handle([&](const path_handle_t& p){
        paths.emplace_back(p);
    });

    unordered_set<handle_t> empty_nodes;

    for (auto& p: paths){
        vector <handle_t> spliced_path;
        bool found_empty_nodes = false;
        string path_name = graph.get_path_name(p);

        graph.for_each_step_in_path(p, [&](const step_handle_t& s){
            auto h = graph.get_handle_of_step(s);

            if (graph.get_length(h) == 0){
                if (empty_nodes.count(h) == 0) {
                    // Short circuit this node by connecting all left neighbors to all right neighbors.
                    // But don't do it if this node was visited before.
                    splice_node_neighbors(graph, h);
                }

                empty_nodes.emplace(h);
                found_empty_nodes = true;
            }
            else{
                // Only add non-empty nodes to the updated path
                spliced_path.emplace_back(h);
            }
        });

        if (found_empty_nodes){
            graph.destroy_path(p);

            // Construct new path that doesn't have any empty nodes
            auto new_path = graph.create_path_handle(path_name);
            for (auto& h: spliced_path){
                graph.append_step(new_path, h);
            }
        }
    }

    graph.for_each_handle([&](const handle_t& h){
        if (graph.get_length(h) == 0) {
            empty_nodes.emplace(h);
        }
    });

    for (auto& h: empty_nodes){
        graph.destroy_handle(h);
    }
}


void clean_gfa(path gfa_path){
    HashGraph graph;
    IncrementalIdMap<string> id_map;
    Overlaps overlaps;

    cerr << "Loading GFA..." << '\n';
    gfa_to_handle_graph(graph, id_map, overlaps, gfa_path);

    splice_out_empty_nodes(graph);

    path output_path = gfa_path;
    output_path = output_path.replace_extension(".clean.gfa");

    cerr << "Writing to GFA: " << output_path << '\n';
    ofstream file(output_path);
    handle_graph_to_gfa(graph, id_map, overlaps, file);
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

    clean_gfa(gfa_path);

    return 0;
}
