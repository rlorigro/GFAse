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

    gfa_to_handle_graph(graph, id_map, gfa_path, false);

    BubbleGraph bubble_graph;

    graph.for_each_handle([&](const handle_t& h0) {
        auto id0 = graph.get_id(h0);

        // If this is an unphased subgraph, check that it is not sharing its phased neighbors with any other
        // subgraphs, by doing a two-edge walk right/left and left/right
        unordered_set<nid_t> second_degree_neighbors;

        // Make sure it isn't possible to visit more than 2 phased bubble sides from this unphased node
        unordered_set<nid_t> left_first_degree_neighbors;
        unordered_set<nid_t> right_first_degree_neighbors;

        graph.follow_edges(h0, true, [&](const handle_t& h1){
            auto id1 = graph.get_id(h1);
            left_first_degree_neighbors.emplace(id1);

            graph.follow_edges(h1, false, [&](const handle_t& h2){
                auto id2 = graph.get_id(h2);

                if (id0 != id2) {
                    second_degree_neighbors.emplace(id2);
                }
            });
        });

        graph.follow_edges(h0, false, [&](const handle_t& h1){
            auto id1 = graph.get_id(h1);
            right_first_degree_neighbors.emplace(id1);

            graph.follow_edges(h1, true, [&](const handle_t& h2){
                auto id2 = graph.get_id(h2);

                if (id0 != id2) {
                    second_degree_neighbors.emplace(id2);
                }
            });
        });

        int32_t result = -1;

        // If there are no second degree neighbors, this unphased subgraph passes
        if (second_degree_neighbors.size() == 1 and right_first_degree_neighbors.size() < 3 and left_first_degree_neighbors.size() < 3){
            auto id_a = graph.get_id(h0);
            auto id_b = graph.get_id(graph.get_handle(*second_degree_neighbors.begin()));

            result = bubble_graph.try_add_bubble(int32_t(id_a), int32_t(id_b));
        }

        cerr << id_map.get_name(graph.get_id(h0)) << ' ' <<  left_first_degree_neighbors.size() << ' ' <<  right_first_degree_neighbors.size() << ' ' << second_degree_neighbors.size() << ' ' << result << '\n';

    });

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
            "Path to GFA containing assembly graph to be phased");

    app.add_option(
            "-o,--output_dir",
            output_dir,
            "Path to (nonexistent) directory where output will be stored")
            ->required();


    CLI11_PARSE(app, argc, argv);

    find_bubbles_in_gfa(output_dir, gfa_path);

    return 0;
}