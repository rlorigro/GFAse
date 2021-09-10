#include "HaplotypePathKmer.hpp"
#include "IncrementalIdMap.hpp"
#include "BinarySequence.hpp"
#include "handle_to_gfa.hpp"
#include "GraphUtility.hpp"
#include "Filesystem.hpp"
#include "CLI11.hpp"

#include <unordered_set>
#include <map>
#include <string>

using gfase::HaplotypePathKmer;
using gfase::IncrementalIdMap;
using gfase::BinarySequence;
using gfase::for_each_connected_component;
using gfase::split_connected_components;
using gfase::handle_graph_to_gfa;
using gfase::print_graph_paths;
using gfase::plot_graph;

using ghc::filesystem::path;

using std::string;
using std::cout;
using std::cerr;
using std::unordered_set;
using std::map;
using std::runtime_error;


void strip_trailing_space(string& line){
    size_t n_trailing_spaces = 0;
    for (auto iter = line.rbegin(); iter != line.rend(); iter++){
        if (isspace(*iter)){
            n_trailing_spaces++;
        }
    }

    line.resize(line.size() - n_trailing_spaces);
}


void test_performance(path kmer_file_path){
    ifstream file(kmer_file_path);
    string line;

    unordered_set <string> kmers;

    while (getline(file, line)){
        if (line[0] == '>'){
            continue;
        }
        else{
            strip_trailing_space(line);
            kmers.emplace(line);
        }
    }
}


int main (int argc, char* argv[]){
    path kmer_file_path;

    CLI::App app{"App description"};

    app.add_option(
            "-i,--input",
            kmer_file_path,
            "Path to FASTA containing kmers to load into memory")
            ->required();

    CLI11_PARSE(app, argc, argv);

    test_performance(kmer_file_path);

    return 0;
}


