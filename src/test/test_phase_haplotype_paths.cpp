#include "Phase.hpp"

using gfase::phase_haplotype_paths;

int main(){
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();

    // Get test VCF path
    path gfa_path = project_directory / "data/simple_chain.gfa";
    path maternal_kmers_path = project_directory / "data/simple_chain_maternal_kmers.fasta";
    path paternal_kmers_path = project_directory / "data/simple_chain_paternal_kmers.fasta";

    phase_haplotype_paths(gfa_path, 6, paternal_kmers_path, maternal_kmers_path, '.');

    return 0;
}
