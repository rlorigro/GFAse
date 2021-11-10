#include "Phase.hpp"

using gfase::phase_k;

int main(){
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();

    // Get test VCF path
    path gfa_path = project_directory / "data/simple_chain.gfa";
    path maternal_kmers_path = project_directory / "data/simple_chain_maternal_kmers.fasta";
    path paternal_kmers_path = project_directory / "data/simple_chain_paternal_kmers.fasta";

    phase_k(gfa_path, 6, paternal_kmers_path, maternal_kmers_path, '.');

    return 0;
}
