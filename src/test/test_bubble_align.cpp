#include "align.hpp"


void infer_bubbles_from_alignment(
        path output_dir,
        path gfa_path,
        double min_similarity,
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

    GfaReader reader(gfa_path);

    // TODO: move this into domain of bdsg graph instead of GFA reader
    string dummy_name = "";
    string dummy_seq = "";
    vector<Sequence> sequences = {Sequence(dummy_name,dummy_seq)};
    reader.for_each_sequence([&](string& name, string& sequence){
        id_map.insert(name);
        sequences.emplace_back(name, sequence);
    });

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

    vector <pair <string,string> > to_be_aligned;

    get_alignment_candidates(
            sequences,
            id_map,
            to_be_aligned,
            output_dir,
            n_threads,
            sample_rate,
            k,
            n_iterations,
            max_hits,
            min_similarity
    );

    ContactGraph alignment_graph;
    ContactGraph symmetrical_alignment_graph;

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
                    ref(sequences),
                    ref(id_map),
                    ref(alignment_graph),
                    min_similarity,
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

    get_best_overlaps(min_similarity, id_map, alignment_graph, symmetrical_alignment_graph);
    write_alignment_results_to_file(id_map, alignment_graph, symmetrical_alignment_graph, output_dir);

    cerr << t << "Done" << '\n';
}


int main (int argc, char* argv[]){
    path gfa_path;
    path output_dir;
    size_t n_threads = 1;
    size_t max_hits = 1;
    double min_similarity;

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
            "-s,--min_similarity",
            min_similarity,
            "Minimum required similarity to consider it in the bubble graph");

    CLI11_PARSE(app, argc, argv);

    infer_bubbles_from_alignment(output_dir, gfa_path, min_similarity, max_hits, n_threads);

    return 0;
}


