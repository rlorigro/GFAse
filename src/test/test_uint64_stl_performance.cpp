#include "FixedBinarySequence.hpp"
#include "Filesystem.hpp"
#include "CLI11.hpp"

#include <unordered_set>
#include <string>

using ghc::filesystem::path;

using std::unordered_set;
using std::ifstream;
using std::string;
using std::cout;
using std::cerr;
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

    unordered_set <uint64_t> hash_table;

    uint64_t i = 0;
    while (getline(file, line)){
        if (line[0] == '>'){
            continue;
        }
        else{
            hash_table.emplace(i);
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


