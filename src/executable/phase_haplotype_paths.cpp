#include "Filesystem.hpp"
#include "Phase.hpp"
#include "CLI11.hpp"

#include "bdsg/hash_graph.hpp"

#include <string>

using ghc::filesystem::path;

using std::string;
using std::cout;
using std::cerr;


int main (int argc, char* argv[]){
    path gfa_path;
    size_t k;
    path paternal_kmers;
    path maternal_kmers;

    CLI::App app{"App description"};

    app.add_option(
            "-i,--input_gfa",
            gfa_path,
            "Path to GFA containing phased non-overlapping segments")
            ->required();

    app.add_option(
            "-k,--kmer_size",
            k,
            "Length of kmer (k) to use")
            ->required();

    app.add_option(
            "-p,--paternal_kmers",
            paternal_kmers,
            "paternal kmers in FASTA format")
            ->required();

    app.add_option(
            "-m,--maternal_kmers",
            maternal_kmers,
            "maternal kmers in FASTA format")
            ->required();

    CLI11_PARSE(app, argc, argv);

    gfase::phase_haplotype_paths(gfa_path, k, paternal_kmers, maternal_kmers);

    return 0;
}
