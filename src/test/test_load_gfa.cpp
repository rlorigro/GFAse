#include "IncrementalIdMap.hpp"
#include "bdsg/hash_graph.hpp"
#include "gfa_to_handle.hpp"
#include "Filesystem.hpp"
#include "CLI11.hpp"

#include <map>
#include <fstream>
#include <string>

using gfase::IncrementalIdMap;
using bdsg::HashGraph;
using ghc::filesystem::path;

using std::ifstream;
using std::string;
using std::cout;
using std::cerr;
using spp::sparse_hash_set;
using std::map;
using std::runtime_error;


void load_gfa(path gfa_path){
    HashGraph graph;
    IncrementalIdMap<string> id_map;

    cerr << "Loading GFA..." << '\n';

    gfa_to_handle_graph(graph, id_map, gfa_path, false);

}


int main (int argc, char* argv[]){
    path kmer_file_path;

    CLI::App app{"App description"};

    app.add_option(
            "-i,--input",
            kmer_file_path,
            "Path to GFA to test loading")
            ->required();

    CLI11_PARSE(app, argc, argv);

    load_gfa(kmer_file_path);

    return 0;
}


