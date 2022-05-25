#include "CLI11.hpp"
#include "PhaseAssign.hpp"

using CLI::App;
using gfase::CigarSummary;
using gfase::assign_phases;


int main(int argc, char* argv[]){
    path output_dir;
    path pat_ref_path;
    path mat_ref_path;
    path query_path;
    string required_prefix;
    size_t n_threads = 1;

    CLI::App app{"App description"};

    app.add_option(
        "--output_dir",
        output_dir,
        "Path to directory which will be created and contain the results of the evaluation. Directory must not exist yet.")
        ->required();

    app.add_option(
        "--pat_ref",
        pat_ref_path,
        "Path to Fasta containing a phased reference.")
        ->required();

    app.add_option(
        "--mat_ref",
        mat_ref_path,
        "Path to Fasta containing a phased reference.")
        ->required();

    app.add_option(
        "--query",
        query_path,
        "Path to Fasta containing phased sequences to align.")
        ->required();

    app.add_option(
        "-p,required_prefix",
        required_prefix,
        "Optionally skip any reads/queries that dont contain this prefix");

    app.add_option(
        "-t,--max_threads",
        n_threads,
        "Maximum number of threads to use");


    CLI11_PARSE(app, argc, argv);

    unordered_map <string, array<CigarSummary,2> > phased_cigar_summaries;
    array <set <string>, 2> phased_contigs;
    map<string,size_t> query_lengths;

    path query_vs_pat_bam;
    path query_vs_mat_bam;

    assign_phases(
        output_dir,
        pat_ref_path,
        mat_ref_path,
        query_vs_pat_bam,
        query_vs_mat_bam,
        query_path,
        required_prefix,
        n_threads,
        phased_cigar_summaries,
        phased_contigs,
        query_lengths,
        true
    );

    return 0;
}
