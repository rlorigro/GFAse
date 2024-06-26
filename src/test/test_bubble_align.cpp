#include "MultiContactGraph.hpp"
#include "gfa_to_handle.hpp"
#include "Sequence.hpp"
#include "Hasher2.hpp"
#include "Timer.hpp"
#include "align.hpp"

using gfase::construct_alignment_graph;
using gfase::gfa_to_handle_graph;
using gfase::MultiContactGraph;
using gfase::HashResult;
using gfase::Sequence;
using gfase::Timer;


void infer_bubbles_from_alignment(
        path output_dir,
        path gfa_path,
        double min_ab_over_a,
        double min_ab_over_b,
        size_t max_hits,
        size_t n_threads
        ){

    Timer t;

    if (exists(output_dir)){
        throw runtime_error("ERROR: output directory exists already");
    }
    else {
        create_directories(output_dir);
    }

    // Id-to-name bimap for reference contigs
    IncrementalIdMap<string> id_map(false);

    HashGraph graph;
    Overlaps overlaps;
    gfa_to_handle_graph(graph, id_map, overlaps, gfa_path);

    // Hashing params
    double sample_rate = 0.04;
    size_t k = 22;
    size_t n_iterations = 6;

    // Only align the top n hits
//    size_t max_hits = 5;

    // Hash results must have at least this percent similarity (A U B)/A, where A is larger.
    // Sequence lengths must be at least this ratio.
    // Resulting alignment coverage must be at least this amount on larger node.
//    double min_similarity = 0.2;

    vector <HashResult> to_be_aligned;

    get_alignment_candidates(
            graph,
            id_map,
            to_be_aligned,
            output_dir,
            n_threads,
            sample_rate,
            k,
            n_iterations,
            max_hits,
            min_ab_over_a,
            min_ab_over_b
    );

    MultiContactGraph alignment_graph;
    MultiContactGraph symmetrical_alignment_graph;

    // Thread-related variables
    atomic<size_t> job_index = 0;
    vector<thread> threads;
    mutex output_mutex;

    // Launch threads
    for (uint64_t n=0; n<n_threads; n++){
        try {
            threads.emplace_back(thread(
                    construct_alignment_graph,
                    ref(to_be_aligned),
                    ref(graph),
                    ref(id_map),
                    ref(alignment_graph),
                    min_ab_over_a,
                    ref(output_mutex),
                    ref(job_index)
            ));
        } catch (const exception &e) {
            cerr << e.what() << "\n";
            exit(1);
        }
    }

    // Wait for threads to finish
    for (auto& n: threads){
        n.join();
    }

    get_best_overlaps(min_ab_over_a, id_map, alignment_graph, symmetrical_alignment_graph);
    write_alignment_results_to_file(id_map, alignment_graph, symmetrical_alignment_graph, output_dir);

    cerr << t << "Done" << '\n';
}


int main (int argc, char* argv[]){
    path gfa_path;
    path output_dir;
    size_t n_threads = 1;
    size_t max_hits = 1;
    double min_ab_over_a;
    double min_ab_over_b;

    CLI::App app{"App description"};

    app.add_option(
            "-g,--gfa",
            gfa_path,
            "Path to GFA containing assembly graph to be phased");

    app.add_option(
            "-o,--output_dir",
            output_dir,
            "Path to (nonexistent) directory where output will be stored")
            ->required();

    app.add_option(
            "-t,--threads",
            n_threads,
            "Maximum number of threads to use");

    app.add_option(
            "-n,--max_hits",
            max_hits,
            "Maximum number of overlap hits to consider");

    app.add_option(
            "-a,--min_ab_over_a",
            min_ab_over_a,
            "Minimum required similarity to consider it in the bubble graph");

    app.add_option(
            "-b,--min_ab_over_b",
            min_ab_over_b,
            "Minimum required similarity to consider it in the bubble graph");

    CLI11_PARSE(app, argc, argv);

    infer_bubbles_from_alignment(output_dir, gfa_path, min_ab_over_a, min_ab_over_b, max_hits, n_threads);

    return 0;
}


