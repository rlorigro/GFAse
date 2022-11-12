#include "Hasher2.hpp"

#include <algorithm>
#include <random>


namespace gfase{


HashResult::HashResult(const string& a, const string& b, double ab_over_a, double ab_over_b):
        a(a),
        b(b),
        ab_over_a(ab_over_a),
        ab_over_b(ab_over_b)
{}


HashResult::HashResult():
        a(),
        b(),
        ab_over_a(0),
        ab_over_b(0)
{}


Hasher2::Hasher2(size_t k, double sample_rate, size_t n_iterations, size_t n_threads):
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


const vector<uint64_t> Hasher2::seeds = {
        2502369103967696952, 7135383713162725540, 1362701453977597027,
        4796741865460292034, 871560157616224954, 4803805340556337874,
        1641952223814618101, 2510544324051818521, 1473066970888840136,
        1577134005447637779, 972008570404626173, 1005686003233212348,
        1240204704125780039, 1557652407674820776, 4544521901418279742,
        1251138144701276783, 8602441838744772705, 32484916947076159,
        1315120142294698015, 2509063378545863014, 1482415985806396390,
        885214842288208064, 1535087442639735647, 6497888167064186043,
        1192965173540970572, 6992204336913461076, 1341474913417665702,
        2625384403220687424, 1699601720342247628, 1332231119554771655,
        664080459378528451, 1057988681310196297, 2272107295286101060,
        6953264339271847999, 1249942556427293808, 8944540437674412670,
        1504916406255515393, 1149201722621576609, 8470764489325659499,
        1195888186965799934, 1451808535441337876, 1533913703589838421,
        5403968531279409728, 1326963194949520018, 370318925754183147,
        963549157246851466, 1640661930493169519, 982070786363093029,
        9033896960103175519, 1367523563246311899, 2559561853484687358,
        2651434670380249965, 1099503378964709489, 514739756668963464,
        1465983026425779258, 3772038222307391067, 1364067351799550025,
        1120347697446236631, 2471601015170514172, 5059946906817196428,
        4640483839834766811, 4956910866326391215, 1538218676201653927,
        9006118601272042222, 1771082815585749522, 1681465739857270927,
        3676604234035115931, 9091606640466583270, 1287150277014289671,
        942177590594404633, 1299657487009447182, 5233542693936021032,
        1173915948459697000, 8759703307818868101, 301590180242423745,
        8073837335570366087, 7136899178665934330, 1414892282437583514,
        1318395012090810332, 1025167099994266395, 1428598750082200673,
        1749293293743715507, 7102823587961144657, 1006232440668911839,
        3932036019385053742, 932163345339352343, 1501418999785722112,
        4202779944578005508, 1715699542033256391, 6103254574080548671,
        1867805817184926607, 1655886200461946312, 4941896307182207749,
        768507597447242884, 4793176833258765644, 3559442299447860782,
        2913348424260394462, 5537057559772751678, 371726285994419264,
        1062935676395772243};


uint64_t Hasher2::hash(const BinarySequence<uint64_t>& kmer, size_t seed_index) const{
    return MurmurHash64A(kmer.sequence.data(), int(kmer.get_byte_length()), seeds[seed_index]);
}


///
/// \param sequence
/// \param i iteration of hashing to compute, corresponding to a hash function
void Hasher2::hash_sequence(const Sequence& sequence, int64_t id, const size_t hash_index) {
    BinarySequence<uint64_t> kmer;

    // Forward iteration
    for (auto& c: sequence.sequence) {
        if (BinarySequence<uint64_t>::base_to_index.at(c) == 4){
            // Reset kmer and don't hash any region with non ACGT chars
            kmer = {};
            continue;
        }

        if (kmer.length < k) {
            kmer.push_back(c);
        } else {
            kmer.shift(c);
            uint64_t h = hash(kmer, hash_index);

            if (h < n_bins){
                auto& m = bin_mutexes.at(h % bin_mutexes.size());
                auto bin_index = h % bins.size();
                auto& bin = bins[bin_index];

                if (bin.size() < max_bin_size + 2){
                    m.lock();
                    bin.emplace(id);
                    m.unlock();
                }
            }
        }
    }

    // Reverse complement iteration
    for (auto iter = sequence.sequence.rbegin(); iter != sequence.sequence.rend(); iter++) {
        if (BinarySequence<uint64_t>::base_to_index.at(*iter) == 4){
            // Reset kmer and don't hash any region with non ACGT chars
            kmer = {};
            continue;
        }

        if (kmer.length < k) {
            kmer.push_back(get_reverse_complement(*iter));
        } else {
            kmer.shift(get_reverse_complement(*iter));
            uint64_t h = hash(kmer, hash_index);

            if (h < n_bins){
                auto& m = bin_mutexes.at(h % bin_mutexes.size());
                auto bin_index = h % bins.size();
                auto& bin = bins[bin_index];

                if (bin.size() < max_bin_size + 2){
                    m.lock();
                    bin.emplace(id);
                    m.unlock();
                }
            }
        }
    }
}


void Hasher2::write_hash_frequency_distribution() const{
    map <size_t, size_t> distribution;

    for (auto& b: bins){
        distribution[b.size()]++;
    }

    for (auto& [size, frequency]: distribution){
        cerr << size << '\t' << frequency << '\n';
    }
}


void Hasher2::hash_sequences(const vector<Sequence>& sequences, atomic<size_t>& job_index, const size_t hash_index){
    size_t i = job_index.fetch_add(1);

    while (i < sequences.size()){
        auto id = id_map.get_id(sequences[i].name);
        hash_sequence(sequences[i], id, hash_index);
        i = job_index.fetch_add(1);
    }
}


void Hasher2::hash(const vector<Sequence>& sequences){
    size_t max_kmers_in_sequence = 0;
    for (auto& sequence: sequences) {
        max_kmers_in_sequence += sequence.size();
    }

    cerr << max_kmers_in_sequence << " possible unique kmers in sequence" << '\n';

    // The maximum observable unique kmers is the total number of kmers in the sequence * the sample rate
    max_kmers_in_sequence = size_t(double(max_kmers_in_sequence) * total_sample_rate);

    cerr << max_kmers_in_sequence << " kmers after downsampling" << '\n';
    cerr << max_kmers_in_sequence * bins_scaling_factor << " bins allocated" << '\n';

    // Build ID map for sequence names
    for(auto& sequence: sequences){
        id_map.try_insert(sequence.name);
    }

    // Aggregate results
    for (size_t h=0; h<n_iterations; h++){
        cerr << "Beginning iteration: " << h << '\n';

        bins.clear();
        bins.resize(max_kmers_in_sequence * bins_scaling_factor);

        for (auto& b: bins){
            bins.reserve(max_bin_size+2);
        }

        // Thread-related variables
        atomic<size_t> job_index = 0;
        vector<thread> threads;

        // Launch threads
        for (uint64_t t=0; t<n_threads; t++){
            try {
                threads.emplace_back(thread(
                        &Hasher2::hash_sequences,
                        this,
                        ref(sequences),
                        ref(job_index),
                        h
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

        // Iterate all hash bins for this iteration (unique hash function)
        for (auto& bin: bins){
            if (bin.size() > max_bin_size){
                continue;
            }

            vector <int64_t> items(bin.size());

            size_t n = 0;
            for (auto& id: bin){
                items[n] = id;
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


void Hasher2::hash(const HandleGraph& graph, const IncrementalIdMap<string>& graph_id_map){
    vector<Sequence> sequences;
    sequences.reserve(graph.get_node_count());

    graph.for_each_handle([&](const handle_t& h){
        auto n = graph.get_id(h);
        auto name = graph_id_map.get_name(n);
        auto sequence = graph.get_sequence(h);
        sequences.emplace_back(name, sequence);
    });

    auto rng = std::default_random_engine(seeds[0]);
    std::shuffle(std::begin(sequences), std::end(sequences), rng);

    hash(sequences);
}


void Hasher2::get_best_matches(map<string, string>& matches, double certainty_threshold) const{
    for (auto& [id, results]: overlaps){
        auto total_hashes = double(results.at(id));

        if (total_hashes < double(min_hashes)){
            continue;
        }

        map <size_t, int64_t> sorted_scores;

        for (auto& [other_id, score]: results){
            // Skip self-hits
            if (other_id == id){
                continue;
            }

            sorted_scores.emplace(score,other_id);
        }

        if (sorted_scores.empty()){
            continue;
        }

        auto max_id = sorted_scores.rbegin()->second;
        auto max_hashes = double(sorted_scores.rbegin()->first);

        // Just take any top hit with greater than % threshold match, later will be used during symmetry filtering
        if (max_hashes/double(total_hashes) > certainty_threshold){
            auto name = id_map.get_name(id);
            auto max_name = id_map.get_name(max_id);
            matches[name] = max_name;
        }
    }
}


void Hasher2::convert_to_contact_graph(
        ContactGraph& contact_graph,
        const IncrementalIdMap<string>& graph_id_map,
        double similarity_threshold,
        size_t minimum_hashes,
        size_t max_overlaps) const {

    for (auto& [id, results]: overlaps){
        auto total_hashes = double(results.at(id));

        // Don't add every result to the graph. Only consider those with at least a certain number of hashes
        if (total_hashes < double(minimum_hashes)){
            continue;
        }

        // For each node try adding it to the graph, and give it a "coverage" that corresponds to its number of hashes
        auto id_a = graph_id_map.get_id(id_map.get_name(id));
        contact_graph.try_insert_node(int32_t(id_a));
        contact_graph.set_node_coverage(int32_t(id_a), int64_t(total_hashes));

        map <size_t, int64_t> sorted_scores;

        for (auto& [other_id, score]: results){
            // Skip self-hits
            if (other_id == id){
                continue;
            }

            sorted_scores.emplace(score,other_id);
        }

        if (sorted_scores.empty()){
            continue;
        }

        size_t i = 0;

        for (auto iter = sorted_scores.rbegin(); iter != sorted_scores.rend(); ++iter){
            auto score = iter->first;
            auto other_id = iter->second;

            // Don't add every result to the graph. Only consider those above a certain similarity
            if (double(score)/double(total_hashes) < similarity_threshold){
                continue;
            }

            auto id_b = graph_id_map.get_id(id_map.get_name(other_id));
            contact_graph.try_insert_node(int32_t(id_b));
            contact_graph.try_insert_edge(int32_t(id_a), int32_t(id_b), int32_t(score));

            i++;

            // Don't let any node add more than x edges to the graph
            if (i == max_overlaps){
                break;
            }
        }
    }
}


void Hasher2::get_symmetrical_matches(map<string, string>& symmetrical_matches, double certainty_threshold) const{
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


int64_t Hasher2::get_intersection_size(const string& a, const string& b) const{
    int64_t intersection = 0;

    auto id_a = id_map.get_id(a);
    auto result_a = overlaps.find(id_a);

    if (result_a != overlaps.end()){
        auto id_b = id_map.get_id(b);
        auto result_b = result_a->second.find(id_b);

        intersection = result_b->second;
    }

    return intersection;
}


void Hasher2::for_each_overlap(
        size_t max_hits,
        double min_similarity,
        const function<void(const string& a, const string& b, int64_t n_hashes, int64_t total_hashes)>& f) const{

    for (auto& [id, results]: overlaps){
        // Self-hit is the total number of hashes the parent sequence had
        auto total_hashes = double(results.at(id));

        if (total_hashes < double(min_hashes)){
            continue;
        }

        map <size_t, int64_t> sorted_scores;

        for (auto& [other_id, score]: results){
            // Skip self-hits
            if (other_id == id){
                continue;
            }

            sorted_scores.emplace(score,other_id);
        }

        if (sorted_scores.empty()){
            continue;
        }

        size_t i = 0;
        // Report the top hits by % similarity for each
        for (auto iter = sorted_scores.rbegin(); iter != sorted_scores.rend(); ++iter){
            auto score = iter->first;
            auto other_id = iter->second;

            float similarity = float(score)/(float(total_hashes) + 1e-12f);

            if (similarity > min_similarity) {
                auto name = id_map.get_name(id);
                auto other_name = id_map.get_name(other_id);
                f(name, other_name, int64_t(score), int64_t(total_hashes));
            }

            i++;
            if (i == max_hits){
                break;
            }
        }
    }
}


void Hasher2::deallocate_bins(){
    bins = {};
}


void Hasher2::write_results(path output_directory) const{
    path overlaps_path = output_directory / "overlaps.csv";

    ofstream overlaps_file(overlaps_path);

    if (not overlaps_file.good() or not overlaps_file.is_open()){
        throw runtime_error("ERROR: could not write file: " + overlaps_path.string());
    }

    overlaps_file << "name" << ',' << "other_name" << ',' << "score" << ',' << "total_hashes" << ',' << "similarity" << '\n';

    for (auto& [id, results]: overlaps){
        size_t total_hashes = results.at(id);

        map <size_t, int64_t> sorted_scores;

        for (auto& [other_id, score]: results){
            // Skip self-hits
            if (other_id == id){
                continue;
            }

            sorted_scores.emplace(score,other_id);
        }

        size_t i = 0;

        // Report the top hits by % Jaccard similarity for each
        for (auto iter = sorted_scores.rbegin(); iter != sorted_scores.rend(); ++iter){
            auto score = iter->first;
            auto other_id = iter->second;

            double similarity = double(score)/double(total_hashes);

            overlaps_file << id_map.get_name(id) << ',' << id_map.get_name(other_id) << ',' << score << ',' << total_hashes << ',' << similarity << '\n';
            i++;

            if (i == 10){
                break;
            }
        }
    }
}

}