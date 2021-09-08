#include "FixedBinarySequence.hpp"
#include "Filesystem.hpp"
#include "sparsepp/spp.h"
#include "CLI11.hpp"

#include <string>

using ghc::filesystem::path;

using std::ifstream;
using std::string;
using std::cout;
using std::cerr;
using spp::sparse_hash_set;
using std::runtime_error;


void test_performance(path kmer_file_path){
    sparse_hash_set <uint8_t> hash_table;

    hash_table.reserve(1000000000);
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


