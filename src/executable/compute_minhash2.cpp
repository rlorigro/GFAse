#include "Hasher2.hpp"
#include "GfaReader.hpp"
#include "Filesystem.hpp"
#include "Sequence.hpp"
#include "CLI11.hpp"


using gfase::Hasher2;
using gfase::Sequence;


int compute_minhash(path gfa_path, path output_directory, double sample_rate, size_t k, size_t n_iterations, size_t n_threads){
    create_directories(output_directory);

    GfaReader reader(gfa_path);
    Hasher2 hasher(k, sample_rate, n_iterations, n_threads);

    vector<Sequence> sequences;
    reader.for_each_sequence([&](string& name, string& sequence){
        sequences.emplace_back(name, sequence);
    });

    hasher.hash(sequences);
    hasher.write_results(output_directory);

    map<string,string> overlaps;

    hasher.get_symmetrical_matches(overlaps, 0.5);

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
