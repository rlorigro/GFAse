#include "Hasher2.hpp"


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


uint64_t Hasher2::hash(const BinarySequence<uint64_t>& kmer, size_t seed_index) const{
    return MurmurHash64A(kmer.sequence.data(), int(kmer.get_byte_length()), seeds[seed_index]);
}


///
/// \param sequence
/// \param i iteration of hashing to compute, corresponding to a hash function
void Hasher2::hash_sequence(const Sequence& sequence, const size_t hash_index) {
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

                m.lock();
                bins[bin_index].emplace(sequence.name);
                m.unlock();
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

                m.lock();
                bins[bin_index].emplace(sequence.name);
                m.unlock();
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
        hash_sequence(sequences[i], hash_index);
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

    // Aggregate results
    for (size_t h=0; h<n_iterations; h++){
        cerr << "Beginning iteration: " << h << '\n';

        bins.clear();
        bins.resize(max_kmers_in_sequence * bins_scaling_factor);

        for (auto& b: bins){
            bins.reserve(max_bin_size);
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


void Hasher2::hash(const HandleGraph& graph, const IncrementalIdMap<string>& id_map){
    vector<Sequence> sequences;
    sequences.reserve(graph.get_node_count());

    graph.for_each_handle([&](const handle_t& h){
        auto n = graph.get_id(h);
        auto name = id_map.get_name(n);
        auto sequence = graph.get_sequence(h);
        sequences.emplace_back(name, sequence);
    });

    hash(sequences);
}


void Hasher2::get_best_matches(map<string, string>& matches, double certainty_threshold) const{
    for (auto& [name, results]: overlaps){
        auto total_hashes = double(results.at(name));

        if (total_hashes < double(min_hashes)){
            continue;
        }

        map <size_t, string> sorted_scores;

        for (auto& [other_name, score]: results){
            // Skip self-hits
            if (other_name == name){
                continue;
            }

            sorted_scores.emplace(score,other_name);
        }

        if (sorted_scores.empty()){
            continue;
        }

        string max_name = sorted_scores.rbegin()->second;
        auto max_hashes = double(sorted_scores.rbegin()->first);

//        double all_hashes = 0;
//
//        size_t i = 0;
//        // Report the top hits by % Jaccard similarity for each
//        for (auto iter = sorted_scores.rbegin(); iter != sorted_scores.rend(); ++iter){
//            auto score = iter->first;
//            auto other_name = iter->second;
//
//            // First item is the maximum score
//            if (i == 0){
//                max_hashes = double(score);
//                max_name = other_name;
//            }
//
//            all_hashes += double(score);
//
//            i++;
//
//            if (i == 10){
//                break;
//            }
//        }

        // Just take any top hit with greater than % threshold match, later will be used during symmetry filtering
        if (max_hashes/double(total_hashes) > certainty_threshold){
            matches[name] = max_name;
        }
    }
}


void Hasher2::convert_to_contact_graph(
        ContactGraph& contact_graph,
        const IncrementalIdMap<string>& id_map,
        double similarity_threshold,
        size_t minimum_hashes,
        size_t max_overlaps) const {

    for (auto& [name, results]: overlaps){
        auto total_hashes = double(results.at(name));

        // Don't add every result to the graph. Only consider those with at least a certain number of hashes
        if (total_hashes < double(minimum_hashes)){
            continue;
        }

        // For each node try adding it to the graph, and give it a "coverage" that corresponds to its number of hashes
        auto id_a = id_map.get_id(name);
        contact_graph.try_insert_node(int32_t(id_a));
        contact_graph.set_node_coverage(int32_t(id_a), int64_t(total_hashes));

        map <size_t, string> sorted_scores;

        for (auto& [other_name, score]: results){
            // Skip self-hits
            if (other_name == name){
                continue;
            }

            sorted_scores.emplace(score,other_name);
        }

        if (sorted_scores.empty()){
            continue;
        }

        size_t i = 0;

        for (auto iter = sorted_scores.rbegin(); iter != sorted_scores.rend(); ++iter){
            auto score = iter->first;
            auto other_name = iter->second;

            // Don't add every result to the graph. Only consider those above a certain similarity
            if (double(score)/double(total_hashes) < similarity_threshold){
                continue;
            }

            auto id_b = id_map.get_id(other_name);
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

    auto result_a = overlaps.find(a);

    if (result_a != overlaps.end()){
        auto result_b = result_a->second.find(b);

        intersection = result_b->second;
    }

    return intersection;
}


void Hasher2::for_each_overlap(
        size_t max_hits,
        double min_similarity,
        const function<void(const string& a, const string& b, int64_t n_hashes, int64_t total_hashes)>& f) const{

    for (auto& [name, results]: overlaps){
        // Self-hit is the total number of hashes the parent sequence had
        auto total_hashes = double(results.at(name));

        if (total_hashes < double(min_hashes)){
            continue;
        }

        map <size_t, string> sorted_scores;

        for (auto& [other_name, score]: results){
            // Skip self-hits
            if (other_name == name){
                continue;
            }

            sorted_scores.emplace(score,other_name);
        }

        if (sorted_scores.empty()){
            continue;
        }

        string max_name = sorted_scores.rbegin()->second;
        auto max_hashes = double(sorted_scores.rbegin()->first);

        size_t i = 0;
        // Report the top hits by % similarity for each
        for (auto iter = sorted_scores.rbegin(); iter != sorted_scores.rend(); ++iter){
            auto score = iter->first;
            auto other_name = iter->second;

            float similarity = float(score)/(float(total_hashes) + 1e-12f);

            if (similarity > min_similarity) {
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

}