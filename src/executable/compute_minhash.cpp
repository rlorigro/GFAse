#include "BinarySequence.hpp"
#include "GfaReader.hpp"
#include "Filesystem.hpp"
#include "Sequence.hpp"
#include "CLI11.hpp"
#include "misc.hpp"
#include "spp.h"

using ghc::filesystem::path;
using ghc::filesystem::exists;
using ghc::filesystem::create_directories;
using gfase::BinarySequence;
using gfase::Sequence;
using gfase::get_reverse_complement;
using spp::sparse_hash_set;
using spp::sparse_hash_map;

#include <unordered_set>
#include <map>
#include <iostream>
#include <ostream>
#include <atomic>
#include <thread>

using std::unordered_set;
using std::map;
using std::numeric_limits;
using std::stringstream;
using std::ostream;
using std::atomic;
using std::thread;
using std::cerr;
using std::min;
using std::max;


// Override the hash function for integer hashing because the hash is provided as the key already
class Hash{
public:
    size_t operator() (uint64_t const& key) const
    {
        return key;
    }
};


class Equal
{
public:
    bool operator() (uint64_t a, uint64_t b) const
    {
        return a == b;
    }
};


// Where to store the names of reads which share hashed-k-mers
using hash_bins_t = sparse_hash_map <uint64_t, unordered_set <string>, Hash, Equal>;

// How to backtrace for each read to its neighbors
using sketches_t = map <string, sparse_hash_set <uint64_t, Hash, Equal> >;

// Ultimately where the results of LSH are stored
using overlaps_t = sparse_hash_map <string, unordered_map <string, int64_t> >;


class HashCluster{
private:
    vector<hash_bins_t> bins_per_iteration;
    vector<sketches_t> sketches_per_iteration;
    overlaps_t overlaps;

    const size_t k;

    size_t n_possible_bins;
    const size_t n_iterations;
    const double total_sample_rate;
    const double iteration_sample_rate;
    const size_t n_threads;

    size_t n_bins;  // Computed from the above values

    static const vector<uint64_t> seeds;

    /// Methods ///
    void hash_sequence(const Sequence& sequence, size_t i);

public:
    HashCluster(size_t k, double sample_rate, size_t n_iterations, size_t n_threads);
    static uint64_t hash(const BinarySequence<uint64_t>& kmer, size_t seed_index);
    void write_hash_frequency_distribution(const hash_bins_t& bins) const;
    void hash_sequences(const vector<Sequence>& sequences, atomic<size_t>& job_index);
    void hash(const vector<Sequence>& sequences);
    void get_best_matches(map<string, string>& matches, double certainty_threshold);
    void get_symmetrical_matches(map<string, string>& symmetrical_matches, double certainty_threshold);
    void write_results(path output_directory);
};


HashCluster::HashCluster(size_t k, double sample_rate, size_t n_iterations, size_t n_threads):
    k(k),
    n_possible_bins(numeric_limits<uint64_t>::max()),
    n_iterations(n_iterations),
    total_sample_rate(sample_rate),
    iteration_sample_rate(sample_rate/double(n_iterations)),
    n_threads(n_threads)
{
    if (k > 32){
        throw runtime_error("ERROR: cannot perform robust 64bit hashing on kmer of length > 32");
    }

    n_bins = round(double(n_possible_bins)*iteration_sample_rate);

    cerr << "Using " << n_bins << " of " << n_possible_bins << " possible bins, for " << n_iterations 
         << " iterations at a rate of " << iteration_sample_rate << '\n';
}


const vector<uint64_t> HashCluster::seeds = {
        2502369103967696952, 7135383713162725540, 13627014539775970274,
        4796741865460292034, 871560157616224954, 4803805340556337874,
        16419522238146181018, 2510544324051818521, 14730669708888401360,
        15771340054476377792, 9720085704046261736, 10056860032332123480,
        12402047041257800396, 15576524076748207768, 4544521901418279742,
        12511381447012767832, 8602441838744772705, 32484916947076159,
        13151201422946980157, 2509063378545863014, 14824159858063963905,
        885214842288208064, 15350874426397356478, 6497888167064186043,
        11929651735409705723, 6992204336913461076, 13414749134176657021,
        2625384403220687424, 1699601720342247628, 13322311195547716555,
        664080459378528451, 10579886813101962970, 2272107295286101060,
        6953264339271847999, 12499425564272938082, 8944540437674412670,
        15049164062555153936, 11492017226215766095, 8470764489325659499,
        1195888186965799934, 1451808535441337876, 15339137035898384211,
        5403968531279409728, 13269631949495200182, 370318925754183147,
        963549157246851466, 16406619304931695195, 9820707863630930290,
        9033896960103175519, 13675235632463118992, 2559561853484687358,
        2651434670380249965, 10995033789647094898, 514739756668963464,
        14659830264257792589, 3772038222307391067, 13640673517995500257,
        11203476974462366311, 2471601015170514172, 5059946906817196428,
        4640483839834766811, 4956910866326391215, 15382186762016539279,
        9006118601272042222, 17710828155857495220, 16814657398572709278,
        3676604234035115931, 9091606640466583270, 12871502770142896716,
        9421775905944046331, 12996574870094471825, 5233542693936021032,
        11739159484596970007, 8759703307818868101, 301590180242423745,
        8073837335570366087, 7136899178665934330, 14148922824375835145,
        1318395012090810332, 10251670999942663955, 14285987500822006731,
        17492932937437155077, 7102823587961144657, 10062324406689118391,
        3932036019385053742, 9321633453393523433, 15014189997857221125,
        4202779944578005508, 1715699542033256391, 6103254574080548671,
        1867805817184926607, 16558862004619463128, 4941896307182207749,
        768507597447242884, 4793176833258765644, 3559442299447860782,
        2913348424260394462, 5537057559772751678, 371726285994419264,
        10629356763957722439};


uint64_t HashCluster::hash(const BinarySequence<uint64_t>& kmer, size_t seed_index){
    return MurmurHash64A(kmer.sequence.data(), int(kmer.get_byte_length()), seeds[seed_index]);
}

///
/// \param sequence
/// \param i iteration of hashing to compute, corresponding to a hash function
void HashCluster::hash_sequence(const Sequence& sequence, size_t i) {
    BinarySequence<uint64_t> kmer;
//    cerr << sequence.name << ' ' << sequence.size();

    // Forward iteration
    for (auto& c: sequence.sequence) {
        if (kmer.length < k) {
            kmer.push_back(c);
        } else {
            kmer.shift(c);
            uint64_t h = hash(kmer, i);

            if (h < n_bins){
                bins_per_iteration[i][h].emplace(sequence.name);
            }
        }
    }

    // Reverse complement iteration
    for (auto iter = sequence.sequence.rbegin(); iter != sequence.sequence.rend(); iter++) {
        if (kmer.length < k) {
            kmer.push_back(get_reverse_complement(*iter));
        } else {
            kmer.shift(get_reverse_complement(*iter));
            uint64_t h = hash(kmer, i);

            if (h < n_bins){
                bins_per_iteration[i][h].emplace(sequence.name);
            }
        }
    }

//    cerr << ' ' << hashes.size() << '\n';
}


void HashCluster::write_hash_frequency_distribution(const hash_bins_t& bins) const{
    map <size_t, size_t> distribution;

    for (auto& [bin_index, bin]: bins){
        distribution[bin.size()]++;
    }

    for (auto& [size, frequency]: distribution){
        cerr << size << '\t' << frequency << '\n';
    }
}


void HashCluster::hash_sequences(const vector<Sequence>& sequences, atomic<size_t>& job_index){
    size_t i=0;
    while (job_index < n_iterations){
        i = job_index.fetch_add(1);

        for (auto& sequence: sequences) {
            hash_sequence(sequence, i);
        }
    }
}


void HashCluster::hash(const vector<Sequence>& sequences){
    bins_per_iteration.resize(n_iterations);
    sketches_per_iteration.resize(n_iterations);

    // Thread-related variables
    atomic<size_t> job_index = 0;
    vector<thread> threads;

    // Launch threads
    for (uint64_t i=0; i<n_threads; i++){
        try {
            threads.emplace_back(thread(
                    &HashCluster::hash_sequences,
                    this,
                    ref(sequences),
                    ref(job_index)
            ));
        } catch (const exception &e) {
            cerr << e.what() << "\n";
            exit(1);
        }
    }

    // Wait for threads to finish
    for (auto& t: threads){
        t.join();
    }

    // Aggregate results
    for (size_t i=0; i<n_iterations; i++){
        auto& bins = bins_per_iteration[i];

        // Iterate all hash bins for this iteration (unique hash function)
        for (auto& [hash,bin]: bins){
            vector <string> items(bin.size());

            size_t n = 0;
            for (auto& name: bin){
                items[n] = name;
                n++;
            }

            // Iterate all combinations of names found in this bin, including self hits, bc they'll be used as a
            // normalization denominator later.
            for (size_t a=0; a<items.size(); a++){
                for (size_t b=a; b<items.size(); b++){
                    overlaps[items[a]][items[b]]++;

                    // Only increment the reciprocal if it's not a self hit
                    if (a != b) {
                        overlaps[items[b]][items[a]]++;
                    }
                }
            }
        }
    }
}


void HashCluster::get_best_matches(map<string, string>& matches, double certainty_threshold){
    for (auto& [name, results]: overlaps){
        size_t total_hashes = results.at(name);

        map <size_t, string> sorted_scores;

        for (auto& [other_name, score]: results){
            // Skip self-hits
            if (other_name == name){
                continue;
            }

            sorted_scores.emplace(score,other_name);
        }

        size_t i = 0;

        string max_name;
        double max_hashes = 0;
        double all_hashes = 0;

        // Report the top hits by % Jaccard similarity for each
        for (auto iter = sorted_scores.rbegin(); iter != sorted_scores.rend(); ++iter){
            auto score = iter->first;
            auto other_name = iter->second;

            // First item is the maximum score
            if (i == 0){
                max_hashes = double(score);
                max_name = other_name;
            }

            all_hashes += double(score);

            i++;

            if (i == 10){
                break;
            }
        }

        // Use a cheap certainty criteria that asks how much of the total matches were part of the max
        if (max_hashes/all_hashes > certainty_threshold){
            matches[name] = max_name;
        }
    }
}


void HashCluster::get_symmetrical_matches(map<string, string>& symmetrical_matches, double certainty_threshold){
    map<string, string> matches;

    get_best_matches(matches, certainty_threshold);

    for (auto& [a,b]: matches){
        auto b_to_a = matches.find(b);

        if (b_to_a != matches.end()){
            if (b_to_a->second == a){
                // Only accept edges which are symmetrical
                symmetrical_matches.emplace(min(a,b), max(a,b));
            }
        }
    }
}


void HashCluster::write_results(path output_directory){
    path overlaps_path = output_directory / "overlaps.csv";

    ofstream overlaps_file(overlaps_path);

    if (not overlaps_file.good() or not overlaps_file.is_open()){
        throw runtime_error("ERROR: could not write file: " + overlaps_path.string());
    }

    overlaps_file << "name" << ',' << "other_name" << ',' << "score" << ',' << "total_hashes" << ',' << "similarity" << '\n';

    for (auto& [name, results]: overlaps){
        size_t total_hashes = results.at(name);

        map <size_t, string> sorted_scores;

        for (auto& [other_name, score]: results){
            // Skip self-hits
            if (other_name == name){
                continue;
            }

            sorted_scores.emplace(score,other_name);
        }

        size_t i = 0;

        // Report the top hits by % Jaccard similarity for each
        for (auto iter = sorted_scores.rbegin(); iter != sorted_scores.rend(); ++iter){
            auto score = iter->first;
            auto other_name = iter->second;

            double similarity = double(score)/double(total_hashes);

            overlaps_file << name << ',' << other_name << ',' << score << ',' << total_hashes << ',' << similarity << '\n';
            i++;

            if (i == 10){
                break;
            }
        }
    }
}


int compute_minhash(path gfa_path, path output_directory, double sample_rate, size_t k, size_t n_iterations, size_t n_threads){
    create_directories(output_directory);

    GfaReader reader(gfa_path);
    HashCluster hasher(k, sample_rate, n_iterations, n_threads);

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
