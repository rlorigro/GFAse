#include "phase_assign.hpp"

using gfase::assign_phases;


int main(){
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();

    path query_path = project_directory / "data/assembly_run11_mismapped_chr20.fasta";
    path mat_ref_path = project_directory / "data/HG002.mat.cur.20211005.S20.fasta";
    path pat_ref_path = project_directory / "data/HG002.pat.cur.20211005.S20.fasta";

    path output_dir = "test";


    array <set <string>, 2> phased_contigs;
    map<string,size_t> query_lengths;

    assign_phases(
            output_dir,
            pat_ref_path,
            mat_ref_path,
            query_path,
            "",
            4,
            phased_contigs,
            query_lengths,
            true
    );

    return 0;
}
