#include "FixedBinarySequence.hpp"
#include "GfaReader.hpp"
#include "Filesystem.hpp"
#include "CLI11.hpp"
#include "spp.h"

using ghc::filesystem::path;
using gfase::FixedBinarySequence;
using spp::sparse_hash_set;
using spp::sparse_hash_map;

#include <unordered_set>
#include <map>
#include <iostream>

using std::unordered_set;
using std::map;
using std::numeric_limits;
using std::stringstream;
using std::cerr;
using std::cerr;


int compute_minhash(path gfa_path){
    GfaReader reader(gfa_path);

    size_t k = 16;
    FixedBinarySequence<uint32_t,1> s;

    double sampling_rate = 0.05;
    size_t n_possible_bins = numeric_limits<uint32_t>::max();
    auto n_bins = size_t(round(double(n_possible_bins) * sampling_rate));

    cerr << "Using " << n_bins << " of " << n_possible_bins << " possible bins, at a rate of " << sampling_rate << '\n';

    vector <unordered_set <string> > bins(n_bins);
    map <string, sparse_hash_set <uint64_t> > sketches;

    reader.for_each_sequence([&](string& name, string& sequence){
        cerr << name << ' ' << sequence.size() << ' ';
        size_t n_found = 0;

        sparse_hash_set<uint64_t> hashes;

        for (auto& c: sequence){
            s.shift(c, k);
            if (s.sequence.front() < n_bins){
                bins.at(s.sequence.front()).emplace(name);
                hashes.emplace(s.sequence.front());
                n_found++;
            }
        }

        cerr << n_found << '\n';
        sketches.emplace(name, hashes);
    });


    for (auto& [name, hashes]: sketches){
        unordered_map <string, int64_t> overlaps;

        if (hashes.size() < 200){
            continue;
        }

        // For each passing hash that this sequence contained
        for (auto& hash: hashes){
            // Iterate all the other sequences that also contained it
            for (auto& hit: bins.at(hash)){
                // Skip self-hits
                if (hit == name){
                    continue;
                }

                overlaps[hit]++;
            }
        }

        for (auto& [other_name, score]: overlaps){
            double similarity = double(score)/double(hashes.size());
            if (similarity > 0.7) {
                cerr << name << '\t' << other_name << '\t' << score << '\t' << hashes.size() << '\t' << similarity << '\n';
            }
        }
    }

    return 0;
}


int main (int argc, char* argv[]){
    path file_path;

    CLI::App app{"App description"};

    app.add_option(
            "-i,--input",
            file_path,
            "Path to GFA")
            ->required();

    CLI11_PARSE(app, argc, argv);

    compute_minhash(file_path);

    return 0;
}



