#include "IncrementalIdMap.hpp"
#include "gfa_to_handle.hpp"
#include "Filesystem.hpp"
#include "CLI11.hpp"

#include "gfakluge.hpp"
#include "bdsg/hash_graph.hpp"

#include <string>

using gfase::IncrementalIdMap;
using gfase::gfa_to_handle_graph;

using ghc::filesystem::path;
using gfak::GFAKluge;
using bdsg::HashGraph;
using bdsg::path_handle_t;
using bdsg::step_handle_t;
using bdsg::handle_t;

using std::string;


void unzip(path gfa_path){
    HashGraph graph;
    IncrementalIdMap<string> id_map;

    gfa_to_path_handle_graph(gfa_path, graph, id_map);

    cout << graph.get_path_count() << '\n';
    graph.for_each_path_handle([&](const path_handle_t& p) {
        cout << graph.get_path_name(p) << "\n";

        graph.for_each_step_in_path(p, [&](const step_handle_t s){
            handle_t h = graph.get_handle_of_step(s);
            nid_t n = graph.get_id(h);
            string sequence = graph.get_sequence(h);
            string name = id_map.get_name(n);

            cout << '\t' << name << " " << sequence << '\n';
        });
    });
}


int main (int argc, char* argv[]){
    path gfa_path;

    CLI::App app{"App description"};

    app.add_option("-i,--input_gfa", gfa_path, "Path to GFA containing phased non-overlapping segments")
            ->required();

    CLI11_PARSE(app, argc, argv);

    unzip(gfa_path);

    return 0;
}
