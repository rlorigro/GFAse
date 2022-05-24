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


void for_each_two_hop_neighbor(){

}


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

        // Do a two-edge walk right/left and left/right
        unordered_set<nid_t> left_second_degree_neighbors;
        unordered_set<nid_t> right_second_degree_neighbors;

        unordered_set<nid_t> left_first_degree_neighbors;
        unordered_set<nid_t> right_first_degree_neighbors;

        graph.follow_edges(h0, true, [&](const handle_t& h1){
            auto id1 = graph.get_id(h1);
            left_first_degree_neighbors.emplace(id1);

            graph.follow_edges(h1, false, [&](const handle_t& h2){
                auto id2 = graph.get_id(h2);

                if (id0 != id2) {
                    left_second_degree_neighbors.emplace(id2);
                }
            });
        });

        graph.follow_edges(h0, false, [&](const handle_t& h1){
            auto id1 = graph.get_id(h1);
            right_first_degree_neighbors.emplace(id1);

            graph.follow_edges(h1, true, [&](const handle_t& h2){
                auto id2 = graph.get_id(h2);

                if (id0 != id2) {
                    right_second_degree_neighbors.emplace(id2);
                }
            });
        });

        cerr << '\n';
        cerr << "left" << '\n';
        for (auto& id: left_first_degree_neighbors){
            cerr << id << ' ' << id_map.get_name(id) << '\n';
        }
        cerr << "right" << '\n';
        for (auto& id: right_first_degree_neighbors){
            cerr << id << ' ' << id_map.get_name(id) << '\n';
        }
        cerr << "left-right neighbors" << '\n';
        for (auto& id: left_second_degree_neighbors){
            cerr << id << ' ' << id_map.get_name(id) << '\n';
        }
        cerr << "right-left neighbors" << '\n';
        for (auto& id: right_second_degree_neighbors){
            cerr << id << ' ' << id_map.get_name(id) << '\n';
        }

        int32_t result = -1;

        bool is_symmetrical_bubble = (right_second_degree_neighbors == left_second_degree_neighbors);
        bool is_diploid_bubble = (right_second_degree_neighbors.size() == 1);
        bool is_chainable = (right_first_degree_neighbors.size() < 3 and left_first_degree_neighbors.size() < 3);

        // If there are no second degree neighbors, this unphased subgraph passes
        if (is_symmetrical_bubble and is_diploid_bubble and is_chainable){
            auto id_a = graph.get_id(h0);
            auto id_b = graph.get_id(graph.get_handle(*left_second_degree_neighbors.begin()));

            try {
                result = bubble_graph.try_add_bubble(int32_t(id_a), int32_t(id_b));
            }
            catch (exception& e){
                cerr << e.what() << '\n';
                cerr << "a: " << id_a << ',' << id_map.get_name(id_b) << '\n';
                cerr << "b: " << id_b << ',' << id_map.get_name(id_b) << '\n';
                throw runtime_error("");
            }
        }

        cerr << id_map.get_name(graph.get_id(h0)) << ' ' <<  left_first_degree_neighbors.size() << ' ' <<  right_first_degree_neighbors.size() << ' ' << left_second_degree_neighbors.size() << ' ' << right_second_degree_neighbors.size() << ' ' << result << '\n';
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