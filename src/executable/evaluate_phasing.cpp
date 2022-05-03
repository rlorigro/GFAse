#include "phase_assign.hpp"
#include "Filesystem.hpp"
#include "CLI11.hpp"

using gfase::evaluate_phasing;
using ghc::filesystem::path;
using CLI::App;


int main (int argc, char* argv[]){
    path output_dir;
    path phase_csv;
    path pat_ref_path;
    path mat_ref_path;
    path query_path;
    string required_prefix;
    size_t n_threads = 1;

    CLI::App app{"App description"};

    app.add_option(
            "-o,--output_dir",
            output_dir,
            "Path to directory which will be created and contain the results of the evaluation. Directory must not exist yet.")
            ->required();

    app.add_option(
            "-c,--phase_csv",
            phase_csv,
            "Path to Fasta containing a phased reference.")
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
            "-q,--query",
            query_path,
            "Path to Fasta containing phased sequences to align.")
            ->required();

    app.add_option(
            "-p,--required_prefix",
            required_prefix,
            "Optionally skip any reads/queries that dont contain this prefix");

    app.add_option(
            "-t,--max_threads",
            n_threads,
            "Maximum number of threads to use");


    CLI11_PARSE(app, argc, argv);

    evaluate_phasing(
            output_dir,
            phase_csv,
            pat_ref_path,
            mat_ref_path,
            query_path,
            required_prefix,
            n_threads,
            true
    );

    return 0;
}