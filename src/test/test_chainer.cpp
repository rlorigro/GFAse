#include "MultiContactGraph.hpp"
#include "IncrementalIdMap.hpp"
#include "gfa_to_handle.hpp"
#include "graph_utility.hpp"
#include "hash_graph.hpp"
#include "Filesystem.hpp"
#include "Chainer.hpp"
#include "CLI11.hpp"

using gfase::gfa_to_handle_graph;
using gfase::MultiContactGraph;
using gfase::IncrementalIdMap;
using gfase::Chainer;

using bdsg::HashGraph;


using ghc::filesystem::path;
using CLI::App;

#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <stdexcept>
#include <fstream>
#include <utility>
#include <limits>
#include <vector>
#include <array>
#include <set>

using std::numeric_limits;
using std::unordered_set;
using std::unordered_map;
using std::runtime_error;
using std::exception;
using std::to_string;
using std::function;
using std::ifstream;
using std::ofstream;
using std::shuffle;
using std::string;
using std::vector;
using std::array;
using std::pair;
using std::stoi;
using std::cerr;
using std::ref;
using std::set;


void chain(path output_dir, path gfa_path){
    if (exists(output_dir)){
        throw runtime_error("ERROR: output directory exists already");
    }
    else {
        create_directories(output_dir);
    }

    // Id-to-name bimap for reference contigs
    IncrementalIdMap<string> id_map(false);
    HashGraph graph;
    MultiContactGraph contact_graph;

    vector<string> phase_0_names = {"b2", "c", "c3", "diploid_a1", "diploid_b1", "diploid_c0", "diploid_d1", "diploid_e1", "e2", "e3", "f", "h", "i2", "i3", "k2_1", "l", "l2_0", "l3", "n", "o2_1", "p2_0", "q", "s", "tip_b", "tip_f", "tree_c0", "tree_c2"};
    vector<string> phase_1_names = {"b", "b3", "c2", "diploid_a0", "diploid_b0", "diploid_c1", "diploid_d0", "diploid_e0", "e", "f2", "f3", "h2", "h3", "i", "k", "k2_0", "k3", "l2_1", "o", "o2_0", "p2_1", "r", "t", "tip_c", "tip_e", "tree_c1", "tree_c3"};

    // Construct graph from GFA
    gfa_to_handle_graph(graph, id_map, gfa_path, false, true);

    for (auto& item: phase_0_names){
        auto id = int32_t(id_map.get_id(item));
        contact_graph.try_insert_node(id);
        contact_graph.set_partition(id, -1);
    }

    for (auto& item: phase_1_names){
        auto id = int32_t(id_map.get_id(item));
        contact_graph.try_insert_node(id);
        contact_graph.set_partition(id, 1);
    }

    Chainer chainer;
    chainer.generate_chain_paths(graph, id_map, contact_graph);

    path chained_gfa_path = output_dir / "chained.gfa";
    ofstream chained_gfa(chained_gfa_path);

    if (not (chained_gfa.is_open() and chained_gfa.good())){
        throw runtime_error("ERROR: could not write to file: " + chained_gfa_path.string());
    }

    handle_graph_to_gfa(graph, id_map, chained_gfa);

    unzip(graph, id_map, false, false);

    chainer.write_chainable_nodes_to_bandage_csv(output_dir, id_map);

    path unzipped_gfa_path = output_dir / "unzipped.gfa";
    ofstream unzipped_gfa(unzipped_gfa_path);

    if (not (unzipped_gfa.is_open() and unzipped_gfa.good())){
        throw runtime_error("ERROR: could not write to file: " + unzipped_gfa_path.string());
    }

    handle_graph_to_gfa(graph, id_map, unzipped_gfa);
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

    chain(output_dir, gfa_path);

    return 0;
}

