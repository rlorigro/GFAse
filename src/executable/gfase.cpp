#include "IncrementalIdMap.hpp"
#include "gfa_to_handle.hpp"
#include "handle_to_gfa.hpp"
#include "GraphUtility.hpp"
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


void unzip(MutablePathDeletableHandleGraph& graph, IncrementalIdMap<string>& id_map){
    unordered_set<nid_t> nodes_to_be_destroyed;
//    unordered_set<path_handle_t> paths_to_be_destroyed;

    cout << graph.get_path_count() << '\n';
    vector<string> path_names;

    vector<path_handle_t> paths;

    graph.for_each_path_handle([&](const path_handle_t& p) {
        paths.emplace_back(p);
    });

    for (auto& p: paths){
        string path_sequence;
        cerr << "Path " << graph.get_path_name(p) << '\n';

        graph.for_each_step_in_path(p, [&](const step_handle_t s){
            handle_t h = graph.get_handle_of_step(s);
            nid_t n = graph.get_id(h);

            string sequence = graph.get_sequence(h);
            path_sequence += sequence;

            string name = id_map.get_name(n);
            cerr << '\t' << name << " " << sequence << '\n';

            nodes_to_be_destroyed.emplace(n);
        });

        // TODO: track provenance, update id_map?
        string haplotype_path_name = graph.get_path_name(p) + "_hap";
        auto new_id = id_map.insert(haplotype_path_name);
        handle_t haplotype_handle = graph.create_handle(path_sequence, new_id);

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

        // Label the new haplotype node using a path to indicate which path it is derived from
        auto haplotype_path_handle = graph.create_path_handle(haplotype_path_name);
        graph.append_step(haplotype_path_handle, haplotype_handle);
    }

//    // Output an image of the graph, can be uncommented for debugging
//    plot_graph(graph, "test_unzip_duplicated");

//    // Destroy all the paths involved in the unzipping
//    for (auto& p: paths){
//        graph.destroy_path(p);
//    }

    // Destroy the nodes that have had their sequences duplicated into haplotypes
    for (auto& n: nodes_to_be_destroyed){
        cerr << "Destroying: " << id_map.get_name(n) << '\n';

//        auto h = graph.get_handle(n);
//        graph.for_each_step_on_handle(h, [&](const step_handle_t& s) {
//            cerr << graph.get_path_name(graph.get_path_handle_of_step(s)) << endl;
//        });

        graph.destroy_handle(graph.get_handle(n));
    }

    // TODO: when creating haplotypes, maintain error bubbles?
    // Search for islands that were created by unphased nodes when the haplotype paths were deleted,
    // and delete the islands (assuming they were errors, and downstream applications won't want bubbles)
    for_each_connected_component(graph, [&](unordered_set<nid_t>& component){
        bool has_path = false;

        // Iterate the connected component and check if it has any haplotype info (if not, it was a unlabeled bubble)
        for (auto& id: component){
            auto h = graph.get_handle(id);

            graph.for_each_step_on_handle(h, [&](const step_handle_t s){
                has_path = true;
                return false;
            });

            cerr << id_map.get_name(graph.get_id(h)) << " has_path: " << has_path << '\n';
        }

        if (not has_path) {
            for (auto& id: component) {
                auto h = graph.get_handle(id);
                graph.destroy_handle(h);
            }
        }
    });
}


void unzip_gfa(path gfa_path){
    HashGraph graph;
    IncrementalIdMap<string> id_map;

    gfa_to_handle_graph(graph, id_map, gfa_path);

    // Output an image of the graph, can be uncommented for debugging
    plot_graph(graph, "test_unzip_unedited");

//    unzip(graph, id_map);

    vector<HashGraph> connected_components;
    vector <IncrementalIdMap<string> > connected_component_ids;

    split_connected_components(graph, id_map, connected_components, connected_component_ids);

    for (size_t i=0; i<connected_components.size(); i++){
        plot_graph(connected_components[i], "component_" + to_string(i));

        cerr << "Component " << to_string(i) << '\n';
        print_graph_paths(connected_components[i], connected_component_ids[i]);

        unzip(connected_components[i], connected_component_ids[i]);
        plot_graph(connected_components[i], "component_" + to_string(i) + "_unzipped");
    }

//    // Output an image of the graph, can be uncommented for debugging
//    plot_graph(graph, "test_unzip_final");
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
