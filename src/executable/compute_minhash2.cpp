#include "Hasher2.hpp"
#include "gfa_to_handle.hpp"
#include "GfaReader.hpp"
#include "Filesystem.hpp"
#include "Sequence.hpp"
#include "CLI11.hpp"


using gfase::gfa_to_handle_graph;
using gfase::HashResult;
using gfase::Hasher2;
using gfase::Sequence;


int compute_minhash(path gfa_path, path output_directory, double sample_rate, size_t k, size_t n_iterations, size_t n_threads){
    create_directories(output_directory);

    HashGraph graph;
    IncrementalIdMap<string> id_map;
    Overlaps overlaps(graph);
    gfa_to_handle_graph(graph, id_map, overlaps, gfa_path);

    Hasher2 hasher(k, sample_rate, n_iterations, n_threads);

    hasher.hash(graph, id_map);
    hasher.write_results(output_directory);

    map<string,string> overlaps;

    unordered_map <pair <string,string>, HashResult> ordered_pairs;

    // Only align the top n hits
    size_t max_hits = 5;

    // Sequence lengths must be at least this ratio.
    // Resulting alignment coverage must be at least this amount on larger node.
    double min_similarity = 0.05;

    // Hash results must have at least this percent similarity (A & B)/A, where A is larger.
    double min_ab_over_a = 0;

    // Hash results must have at least this percent similarity (A & B)/B, where A is larger.
    double min_ab_over_b = 0.7;

    hasher.for_each_overlap(max_hits, min_ab_over_a,[&](const string& a, const string& b, int64_t n_hashes, int64_t total_hashes){
        // Skip self hits
        if (a == b){
            return;
        }

        auto seq_a = graph.get_sequence(graph.get_handle(id_map.get_id(a)));
        auto seq_b = graph.get_sequence(graph.get_handle(id_map.get_id(b)));

        pair<string,string> ordered_pair;
        if (seq_a.size() > seq_b.size()){
            ordered_pair = {a,b};

            auto& result = ordered_pairs[ordered_pair];
            result.ab_over_a = double(n_hashes)/double(total_hashes);
        }
        else{
            ordered_pair = {b,a};

            auto& result = ordered_pairs[ordered_pair];
            result.ab_over_b = double(n_hashes)/double(total_hashes);
        }
    });

    // Filter once both directional hash similarities are established
    for (const auto& [edge,result]: ordered_pairs){
        if (result.ab_over_a >= min_ab_over_a and result.ab_over_b >= min_ab_over_b) {
            cerr << edge.first << ',' << edge.second << ',' << result.ab_over_a << ',' << min_ab_over_a << ',' << result.ab_over_b << ',' << min_ab_over_b << '\n';
        }
    }

    path output_path = output_directory / "pairs.csv";
    ofstream file(output_path);

    if (not file.good() or not file.is_open()){
        throw runtime_error("ERROR: could not write file: " + output_path.string());
    }

    file << "Name" << ',' << "Match" << ',' << "Color" << '\n';
    for (auto& [a,b]: overlaps){
        file << a << ',' << b << ',' << "Cornflower Blue" << '\n';
        file << b << ',' << a << ',' << "Tomato" << '\n';
    }

    return 0;
}


int main (int argc, char* argv[]){
    path file_path;
    path output_directory;
    double sample_rate = 0.1;
    size_t k = 22;
    size_t n_iterations = 10;
    size_t n_threads = 1;

    CLI::App app{"App description"};

    app.add_option(
            "-g,--gfa",
            file_path,
            "Path to GFA")
            ->required();

    app.add_option(
            "-o,--output_directory",
            output_directory,
            "Path to directory where output should be written")
            ->required();

    app.add_option(
            "-r,--sample_rate",
            sample_rate,
            "Sample rate. Proportion [0-1] of k-mers to retain during in comparison")
            ->required();

    app.add_option(
            "-k,--kmer_length",
            k,
            "Length of k-mer to use for hashing")
            ->required();

    app.add_option(
            "-i,--n_iterations",
            n_iterations,
            "Number of iterations (different hash functions), each at a rate of sample_rate/n_iterations")
            ->required();

    app.add_option(
            "-t,--n_threads",
            n_threads,
            "Maximum number of threads to use")
            ->required();

    CLI11_PARSE(app, argc, argv);

    if (sample_rate < 0 or sample_rate > 1){
        throw std::runtime_error("ERROR: sample rate must be between 0 and 1.0");
    }

    compute_minhash(file_path, output_directory, sample_rate, k, n_iterations, n_threads);

    return 0;
}
