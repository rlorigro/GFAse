#include "PhaseAssign.hpp"
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
    path pat_bam_path;
    path mat_bam_path;
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
            "Path to Fasta containing a phased paternal reference.");

    app.add_option(
            "--mat_ref",
            mat_ref_path,
            "Path to Fasta containing a phased maternal reference");

    app.add_option(
            "--pat_bam",
            pat_bam_path,
            "Path to alignment file containing all contigs in assembly aligned to a phased paternal reference. MUST be aligned with explicit mismatch cigars (= and X operations).");

    app.add_option(
            "--mat_bam",
            mat_bam_path,
            "Path to alignment file containing all contigs in assembly aligned to a phased maternal reference. MUST be aligned with explicit mismatch cigars (= and X operations).");

    app.add_option(
            "-q,--query",
            query_path,
            "Path to Fasta containing phased assembly sequences to align.")
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

    bool has_bams = (not pat_bam_path.empty()) and (not mat_bam_path.empty());
    bool has_refs = (not pat_ref_path.empty()) and (not mat_ref_path.empty());

    bool valid_bam_args = has_bams and (not has_refs);
    bool valid_ref_args = (not has_bams) and has_refs;

    if (not (valid_bam_args or valid_ref_args)){
        throw runtime_error("ERROR: commandline arguments must provide a pair of reference sequences (FASTA) or a pair of alignments (BAM) but not both");
    }

    evaluate_phasing(
            output_dir,
            phase_csv,
            pat_ref_path,
            mat_ref_path,
            pat_bam_path,
            mat_bam_path,
            query_path,
            required_prefix,
            n_threads,
            true
    );

    return 0;
}