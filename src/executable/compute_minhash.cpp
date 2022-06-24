#include "BinarySequence.hpp"
#include "GfaReader.hpp"
#include "Filesystem.hpp"
#include "CLI11.hpp"
#include "spp.h"

using ghc::filesystem::path;
using gfase::BinarySequence;
using spp::sparse_hash_set;
using spp::sparse_hash_map;

#include <unordered_set>
#include <map>
#include <iostream>
#include <ostream>

using std::unordered_set;
using std::map;
using std::numeric_limits;
using std::stringstream;
using std::ostream;
using std::cerr;


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


class HashCluster{
private:
    // Where to store the names od reads which share hashed-k-mers
    sparse_hash_map <uint64_t, unordered_set <string>, Hash, Equal> bins;

    // How to backtrace for each read to its neighbors
    map <string, sparse_hash_set <uint64_t, Hash, Equal> > sketches;

    // Ultimately where the results are stored
    sparse_hash_map <string, unordered_map <string, int64_t> > overlaps;

    size_t k;

    double total_sample_rate;
    double iteration_sample_rate;
    size_t n_possible_bins;
    size_t n_iterations;
    size_t n_bins;

    static const vector<uint64_t> seeds;

    /// Methods ///
    void add_sequence(const string& name, const string& sequence, size_t iteration_index);

public:
    HashCluster(size_t k, double sample_rate, size_t n_iterations);
    static uint64_t hash(const BinarySequence<uint64_t>& kmer, size_t seed_index);
    void write_hash_frequency_distribution() const;
    void cluster(GfaReader& reader);
    void write_results();
};


HashCluster::HashCluster(size_t k, double sample_rate, size_t n_iterations):
    k(k),
    total_sample_rate(sample_rate),
    n_iterations(n_iterations)
{
    if (k > 32){
        throw runtime_error("ERROR: cannot perform robust 64bit hashing on kmer of length > 32");
    }
    n_possible_bins = numeric_limits<uint64_t>::max();
    iteration_sample_rate = total_sample_rate/double(n_iterations);
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


void HashCluster::add_sequence(const string& name, const string& sequence, size_t iteration_index) {
    BinarySequence<uint64_t> kmer;
    sparse_hash_set<uint64_t, Hash, Equal> hashes;

    cerr << name << ' ' << sequence.size();

    for (auto& c: sequence) {
        if (kmer.length < k) {
            kmer.push_back(c);
        } else {
            kmer.shift(c);
            uint64_t h = hash(kmer, iteration_index);

            if (h < n_bins){
                bins[h].emplace(name);
                hashes.emplace(h);
            }
        }
    }

    cerr << ' ' << hashes.size() << '\n';

    sketches.emplace(name, hashes);
}


void HashCluster::write_hash_frequency_distribution() const{
    map <size_t, size_t> distribution;

    for (auto& [bin_index, bin]: bins){
        distribution[bin.size()]++;
    }

    for (auto& [size, frequency]: distribution){
        cerr << size << '\t' << frequency << '\n';
    }
}


void HashCluster::cluster(GfaReader& reader){

    for (size_t i=0; i<n_iterations; i++){
        sketches.clear();
        bins.clear();

        reader.for_each_sequence([&](string& name, string& sequence){
            add_sequence(name, sequence, i);
        });

        write_hash_frequency_distribution();

        for (auto& [name, hashes]: sketches){
            // TODO: automate this parameter selection ?
            if (hashes.size() < 10){
                continue;
            }

            // For each passing hash that this sequence contained
            for (auto& hash: hashes){
                // Iterate all the other sequences that also contained it
                for (auto& hit: bins.at(hash)){
                    overlaps[name][hit]++;
                }
            }
        }
    }
}


void HashCluster::write_results(){
    // TODO: sort and filter more comprehensively

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
        for (auto iter = sorted_scores.rbegin(); iter != sorted_scores.rend(); ++iter){
            auto score = iter->first;
            auto other_name = iter->second;

            double similarity = double(score)/double(total_hashes);

            cerr << name << '\t' << other_name << '\t' << score << '\t' << total_hashes << '\t' << similarity << '\n';
            i++;

            if (i == 10){
                break;
            }
        }
    }
}


int compute_minhash(path gfa_path){
    size_t k = 22;

    GfaReader reader(gfa_path);
    HashCluster clusterer(k, 0.12, 3);

    clusterer.cluster(reader);
    clusterer.write_results();

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
