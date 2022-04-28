#include "bamtools/api/BamReader.h"
#include "IncrementalIdMap.hpp"
#include "Filesystem.hpp"
#include "sparsepp/spp.h"
#include "CLI11.hpp"
#include "misc.hpp"
#include "Sam.hpp"

using gfase::for_element_in_sam_file;
using gfase::IncrementalIdMap;
using gfase::SamElement;

using BamTools::BamAlignment;
using BamTools::BamReader;

using ghc::filesystem::path;
using spp::sparse_hash_map;
using CLI::App;

#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <stdexcept>
#include <algorithm>
#include <fstream>
#include <utility>
#include <atomic>
#include <thread>
#include <random>
#include <limits>
#include <bitset>
#include <vector>
#include <mutex>
#include <array>
#include <set>

using std::numeric_limits;
using std::unordered_set;
using std::unordered_map;
using std::runtime_error;
using std::streamsize;
using std::exception;
using std::to_string;
using std::function;
using std::ifstream;
using std::ofstream;
using std::shuffle;
using std::string;
using std::vector;
using std::bitset;
using std::thread;
using std::atomic;
using std::array;
using std::mutex;
using std::pair;
using std::stoi;
using std::cerr;
using std::cref;
using std::ref;
using std::set;


using paired_mappings_t = sparse_hash_map <string, array <set <SamElement>, 2> >;
using unpaired_mappings_t = sparse_hash_map <string, set <SamElement> >;
using contact_map_t = sparse_hash_map <int32_t, sparse_hash_map<int32_t, int32_t> >;


void print_mappings(const paired_mappings_t& mappings){
    for (const auto& [name,mates]: mappings){
        cerr << '\n';
        cerr << name << '\n';

        cerr << "First mates:" << '\n';
        for (auto& e: mates[0]){
            cerr << '\t' << e.ref_name << '\t' << bitset<sizeof(e.flag)*8>(e.flag) << '\n';
            cerr << '\t' << e.is_first_mate() << ' ' << e.is_second_mate() << '\n';
        }

        cerr << "Second mates:" << '\n';
        for (auto& e: mates[1]){
            cerr << '\t' << e.ref_name << '\t' << bitset<sizeof(e.flag)*8>(e.flag) << '\n';
            cerr << '\t' << e.is_first_mate() << ' ' << e.is_second_mate() << '\n';
        }
    }
}


void print_mappings(const unpaired_mappings_t& mappings){
    for (const auto& [name,elements]: mappings){
        cerr << '\n';
        cerr << name << '\n';

        cerr << "Alignments:" << '\n';
        for (auto& e: elements){
            cerr << '\t' << e.ref_name << '\t' << bitset<sizeof(e.flag)*8>(e.flag) << '\t' << "MQ" << int(e.mapq) << ' ' << 'P' << !e.is_not_primary() << ' ' << 'S' << e.is_supplementary() << '\n';
        }
    }
}


class Bubble{
public:
    array<int32_t,2> ids;
    bool phase;

    Bubble();
    Bubble(int32_t id1, int32_t id2, bool phase);
    void flip();
    int32_t first() const;
    int32_t second() const;
    int32_t get(bool side) const;
    int32_t is_first(int32_t id) const;
    int32_t is_second(int32_t id) const;
};


Bubble::Bubble(int32_t id1, int32_t id2, bool phase):
        ids({id1,id2}),
        phase(phase)
{}


Bubble::Bubble():
        ids({-1,-1}),
        phase(0)
{}


void Bubble::flip(){
    phase = not phase;
}


int32_t Bubble::first() const{
    return ids[0 + phase];
}


int32_t Bubble::second() const{
    return ids[1 - phase];
}

///
/// s p return
/// 0 1 1
/// 1 0 1
/// 0 0 0
/// 1 1 0
int32_t Bubble::get(bool side) const{
    return ids[side != phase];
}


int32_t Bubble::is_first(int32_t id) const{
    if (get(id) == 0){
        return true;
    }
    else{
        return false;
    }
}


int32_t Bubble::is_second(int32_t id) const{
    if (get(id) == 1){
        return true;
    }
    else{
        return false;
    }
}


class Bubbles {
public:
    vector<Bubble> bubbles;
    unordered_map<int32_t,int32_t> ref_id_to_bubble_id;
    vector <vector <int32_t> > bubble_to_bubble;
    unordered_set <pair <int32_t, int32_t> > bubble_pairs;

    Bubbles();
    void write_bandage_csv(path output_path, IncrementalIdMap<string>& id_map);
    void generate_bubble_adjacency_from_contact_map(const contact_map_t& contact_map);
    void for_each_adjacent_bubble(int32_t b, const function<void(Bubble& bubble)>& f);
    void for_each_adjacent_bubble(int32_t b, const function<void(const Bubble& bubble)>& f) const;
    void for_each_bubble_pair(const function<void(Bubble& b1, Bubble& b2)>& f);
    void for_each_bubble_pair(const function<void(const Bubble& b1, const Bubble& b2)>& f) const;
    void get_phases(vector<bool>& bubble_phases) const;
    void set_phases(const vector<bool>& bubble_phases);
    void at(size_t i, Bubble& b);
    Bubble at(size_t i) const;
    void emplace(int32_t id1, int32_t id2, bool phase);
    int32_t find(int32_t id);
    void flip(size_t b);
    size_t size() const;
};


Bubbles::Bubbles():
        bubbles(),
        ref_id_to_bubble_id(),
        bubble_to_bubble(),
        bubble_pairs()
{}


void Bubbles::generate_bubble_adjacency_from_contact_map(const contact_map_t& contact_map){
    bubble_to_bubble.resize(bubbles.size());

    for (size_t b=0; b<bubbles.size(); b++) {
        auto id0 = bubbles[b].get(0);
        auto id1 = bubbles[b].get(1);

        unordered_set<size_t> other_bubbles;

        auto c0 = contact_map.find(id0);
        auto c1 = contact_map.find(id1);

        // Iterate all contacts with id0
        if (c0 != contact_map.end()) {
            for (auto&[other_id, count]: c0->second) {
                auto b_other = find(other_id);
                other_bubbles.emplace(b_other);
            }
        }

        // Iterate all contacts with id1
        if (c1 != contact_map.end()) {
            for (auto&[other_id, count]: c1->second) {
                auto b_other = find(other_id);
                other_bubbles.emplace(b_other);
            }
        }

        for (auto& b_other: other_bubbles) {
            bubble_to_bubble[b].emplace_back(b_other);
            bubble_pairs.emplace(b,b_other);
        }
    }
}


void Bubbles::for_each_adjacent_bubble(int32_t b, const function<void(Bubble& bubble)>& f){
    for (auto& b_other: bubble_to_bubble[b]){
        f(bubbles[b_other]);
    }
}


void Bubbles::for_each_adjacent_bubble(int32_t b, const function<void(const Bubble& bubble)>& f) const{
    for (const auto& b_other: bubble_to_bubble[b]){
        f(bubbles[b_other]);
    }
}


void Bubbles::for_each_bubble_pair(const function<void(Bubble& b0, Bubble& b1)>& f){
    for (auto& p: bubble_pairs){
        f(bubbles[p.first], bubbles[p.second]);
    }
}


void Bubbles::for_each_bubble_pair(const function<void(const Bubble& b0, const Bubble& b1)>& f) const{
    for (const auto& p: bubble_pairs){
        f(bubbles[p.first], bubbles[p.second]);
    }
}


void Bubbles::write_bandage_csv(path output_path, IncrementalIdMap<string>& id_map){
    ofstream file(output_path);

    if (not file.is_open() or not file.good()){
        throw std::runtime_error("ERROR: could not write to file: " + output_path.string());
    }

    file << "Name" << ',' << "Phase" << ',' << "Color" << '\n';

    for (auto& b: bubbles){
        auto id0 = b.get(0);
        auto id1 = b.get(1);

        file << id_map.get_name(id0) << ',' << 0 << ',' << "Dark Orange" << '\n';
        file << id_map.get_name(id1) << ',' << 1 << ',' << "Green Yellow" << '\n';
    }
}


void Bubbles::emplace(int32_t id1, int32_t id2, bool phase){
    ref_id_to_bubble_id[id1] = int32_t(bubbles.size());
    ref_id_to_bubble_id[id2] = int32_t(bubbles.size());
    bubbles.emplace_back(id1, id2, phase);
}


void Bubbles::get_phases(vector<bool>& bubble_phases) const{
    if (bubble_phases.size() != bubbles.size()){
        bubble_phases.resize(bubbles.size());
    }

    for (size_t b=0; b<bubbles.size(); b++){
        bubble_phases[b] = bubbles[b].phase;
    }
}


void Bubbles::set_phases(const vector<bool>& bubble_phases){
    if (bubble_phases.size() != bubbles.size()){
        throw runtime_error("ERROR: cannot set phase vector with unequal size vector");
    }

    for (size_t b=0; b<bubbles.size(); b++){
        bubbles[b].phase = bubble_phases[b];
    }
}


size_t Bubbles::size() const{
    return bubbles.size();
}


void Bubbles::at(size_t i, Bubble& b){
    b = bubbles.at(i);
}


Bubble Bubbles::at(size_t i) const{
    return bubbles.at(i);
}


int32_t Bubbles::find(int32_t id){
    return ref_id_to_bubble_id.at(id);
}


void Bubbles::flip(size_t b){
    bubbles.at(b).flip();
}


void generate_bubbles_from_shasta_names(Bubbles& bubbles, IncrementalIdMap<string>& id_map){
    unordered_set<int32_t> visited;

    for (auto& [name,id]: id_map.ids){
        if (visited.count(int32_t(id)) > 0){
            continue;
        }

        // Split to find last field, which should be 0/1 for shasta PR segments
        auto i = name.find_last_of('.');

        if (i < 2){
            throw std::runtime_error("ERROR: shasta phase field unconventional: " + name);
        }

        auto prefix = name.substr(0, i);
        int64_t side = stoi(name.substr(i+1));

        // Cheap test to check for proper syntax
        if (side > 1 or side < 0){
            throw std::runtime_error("ERROR: shasta bubble side not 0/1: " + name);
        }

        // Find complement
        string other_name = prefix + '.' + to_string(1 - side);

        // Look for the other name in the id_map, add it if it doesn't exist
        // Later, will need to assume these missing entries in the contact map are 0
        auto other_id = id_map.try_insert(other_name);
        bubbles.emplace(int32_t(id), int32_t(other_id), 0);

        visited.emplace(id);
        visited.emplace(other_id);
    }
}


void parse_paired_sam_file(
        path sam_path,
        paired_mappings_t& mappings,
        IncrementalIdMap<string>& id_map,
        string required_prefix,
        int8_t min_mapq){

    for_element_in_sam_file(sam_path, [&](SamElement& e){
        bool valid_prefix = true;
        if (not required_prefix.empty()){
            for (size_t i=0; i<required_prefix.size(); i++){
                if (e.ref_name[i] != required_prefix[i]){
                    valid_prefix = false;
                    break;
                }
            }
        }

        if (valid_prefix) {
            id_map.try_insert(e.ref_name);

            if (e.mapq >= min_mapq and (not e.is_not_primary())) {
                mappings[e.read_name][e.is_second_mate()].emplace(e);
            }
        }
    });
}


void parse_unpaired_sam_file(
        path sam_path,
        unpaired_mappings_t& mappings,
        IncrementalIdMap<string>& id_map,
        string required_prefix,
        int8_t min_mapq){

    for_element_in_sam_file(sam_path, [&](SamElement& e){
        bool valid_prefix = true;
        if (not required_prefix.empty()){
            for (size_t i=0; i<required_prefix.size(); i++){
                if (e.ref_name[i] != required_prefix[i]){
                    valid_prefix = false;
                    break;
                }
            }
        }

        if (valid_prefix) {
            id_map.try_insert(e.ref_name);

            if (e.mapq >= min_mapq and (not e.is_not_primary())) {
                mappings[e.read_name].emplace(e);
            }
        }
    });
}


void parse_unpaired_bam_file(
        path bam_path,
        unpaired_mappings_t& mappings,
        IncrementalIdMap<string>& id_map,
        string required_prefix,
        int8_t min_mapq){

    BamReader reader;

    if (!reader.Open(bam_path) ) {
        throw std::runtime_error("ERROR: could not read BAM file: " + bam_path.string());
    }

    auto reference_data = reader.GetReferenceData();

    BamAlignment e;
    size_t l = 0;

    while (reader.GetNextAlignment(e) ) {
//        cerr << "e.Name" << ' ' << e.Name << '\n';
//        cerr << "e.RefID" << ' ' << e.RefID << '\n';
//        cerr << "int32_t(l)" << ' ' << int32_t(l) << '\n';
//        cerr << "int16_t(e.AlignmentFlag)" << ' ' << int(e.AlignmentFlag) << '\n';
//        cerr << "int8_t(e.MapQuality)" << ' ' << int(e.MapQuality) << '\n';
//        cerr << '\n';

        auto& ref_name = reference_data.at(e.RefID).RefName;

        bool valid_prefix = true;
        if (not required_prefix.empty()){
            for (size_t i=0; i<required_prefix.size(); i++){
                if (ref_name[i] != required_prefix[i]){
                    valid_prefix = false;
                    break;
                }
            }
        }

        if (valid_prefix) {
            id_map.try_insert(ref_name);

            if (e.MapQuality >= min_mapq and e.IsPrimaryAlignment()) {
                SamElement s(e.Name, ref_name, int32_t(l), int16_t(e.AlignmentFlag), int8_t(e.MapQuality));
                mappings[e.Name].emplace(s);
            }
        }

        l++;
    }
}


void remove_unpaired_reads(paired_mappings_t& mappings){
    vector<string> to_be_deleted;
    for (auto& [name,mates]: mappings){

        bool has_first_mate = not mates[0].empty();
        bool has_second_mate = not mates[1].empty();

        if (not (has_first_mate and has_second_mate)){
            to_be_deleted.emplace_back(name);
        }
    }

    for (auto& item: to_be_deleted){
        mappings.erase(item);
    }
}


void remove_singleton_reads(unpaired_mappings_t& mappings){
    vector<string> to_be_deleted;
    for (auto& [name,elements]: mappings){

        // Remove all the entries where only one mapping exists for that read
        if (elements.size() < 2){
            to_be_deleted.emplace_back(name);
        }
    }

    for (auto& item: to_be_deleted){
        mappings.erase(item);
    }
}


void generate_contact_map_from_mappings(
        unpaired_mappings_t& mappings,
        IncrementalIdMap<string>& id_map,
        contact_map_t& contact_map
){
    // All-vs-all, assuming reference names are unique in first vs second mates
    for (auto& [name,elements]: mappings) {
        size_t i = 0;

        for (auto& e: elements){
            size_t j = 0;

            for (auto& e2: elements){
                if (j > i){
                    if (e.ref_name != e2.ref_name) {
                        auto id = int32_t(id_map.get_id(e.ref_name));
                        auto id2 = int32_t(id_map.get_id(e2.ref_name));

                        contact_map[id][id2]++;
                        contact_map[id2][id]++;
                    }
                }
                j++;
            }
            i++;
        }
    }
}


void generate_contact_map_from_mappings(
        paired_mappings_t& mappings,
        IncrementalIdMap<string>& id_map,
        contact_map_t& contact_map
){
    // All-vs-all, assuming reference names are unique in first vs second mates
    for (auto& [name,mates]: mappings) {
        for (auto& e: mates[0]){
            for (auto& e2: mates[1]){
                if (e.ref_name != e2.ref_name) {
                    auto id = int32_t(id_map.get_id(e.ref_name));
                    auto id2 = int32_t(id_map.get_id(e2.ref_name));

                    contact_map[id][id2]++;
                    contact_map[id2][id]++;
                }
            }
        }
    }
}


void generate_adjacency_matrix(
        const Bubbles& bubbles,
        const contact_map_t& contact_map,
        vector <vector <int32_t> >& adjacency
){

    size_t n = 2*bubbles.size();
    adjacency.resize(n, vector<int32_t>(n, 0));

    for (auto& [id1,map2]: contact_map){
        for (auto& [id2,count]: map2){
            adjacency[id1][id2] = count;
        }
    }
}


int64_t compute_total_consistency_score(
        const Bubbles& bubbles,
        const contact_map_t& contact_map
){
    int64_t score = 0;

    bubbles.for_each_bubble_pair([&](const Bubble& b0, const Bubble& b1){
        auto id0 = b0.get(0);
        auto id1 = b0.get(1);

        auto other_id0 = b1.get(0);
        auto other_id1 = b1.get(1);

        auto r0 = contact_map.find(id0);
        if (r0 != contact_map.end()){
            auto r00 = r0->second.find(other_id0);
            auto r01 = r0->second.find(other_id1);

            if (r00 != r0->second.end()){
                score += r00->second;
            }

            if (r01 != r0->second.end()){
                score -= r01->second;
            }
        }

        auto r1 = contact_map.find(id1);
        if (r1 != contact_map.end()){
            auto r10 = r1->second.find(other_id0);
            auto r11 = r1->second.find(other_id1);

            if (r10 != r1->second.end()){
                score -= r10->second;
            }

            if (r11 != r1->second.end()){
                score += r11->second;
            }
        }
    });

    return score;
}


int64_t compute_consistency_score(
        const Bubbles& bubbles,
        size_t bubble_index,
        const contact_map_t& contact_map
        ){
    int64_t score = 0;

    const Bubble bubble = bubbles.at(int32_t(bubble_index));

    auto id0 = bubble.get(0);
    auto id1 = bubble.get(1);

    bubbles.for_each_adjacent_bubble(int32_t(bubble_index), [&](const Bubble& other_bubble){
        auto other_id0 = other_bubble.get(0);
        auto other_id1 = other_bubble.get(1);

        auto r0 = contact_map.find(id0);
        if (r0 != contact_map.end()){
            auto r00 = r0->second.find(other_id0);
            auto r01 = r0->second.find(other_id1);

            if (r00 != r0->second.end()){
                score += r00->second;
            }

            if (r01 != r0->second.end()){
                score -= r01->second;
            }
        }

        auto r1 = contact_map.find(id1);
        if (r1 != contact_map.end()){
            auto r10 = r1->second.find(other_id0);
            auto r11 = r1->second.find(other_id1);

            if (r10 != r1->second.end()){
                score -= r10->second;
            }

            if (r11 != r1->second.end()){
                score += r11->second;
            }
        }
    });

    return score;
}


void random_phase_search(
        const contact_map_t& contact_map,
        const IncrementalIdMap<string>& id_map,
        Bubbles bubbles,
        vector<bool>& best_phases,
        atomic<int64_t>& best_score,
        atomic<size_t>& job_index,
        mutex& phase_mutex,
        size_t m_iterations
        ){

    size_t m = job_index.fetch_add(1);

    // True random number
    std::random_device rd;

    // Pseudorandom generator with true random seed
    std::mt19937 rng(rd());
    std::uniform_int_distribution<int> uniform_distribution(0,int(bubbles.size()-1));

    int64_t total_score;

    vector<size_t> order(bubbles.size());
    for (size_t i=0; i<order.size(); i++){
        order[i] = i;
    }

    while (m < m_iterations) {
        // Randomly perturb
        for (size_t i=0; i < ((bubbles.size()/10) + 1); i++) {
            bubbles.flip(uniform_distribution(rng));
        }

        for (size_t i=0; i < bubbles.size()*3; i++) {
            auto b = uniform_distribution(rng);

            auto score = compute_consistency_score(bubbles, b, contact_map);
            bubbles.flip(b);
            auto flipped_score = compute_consistency_score(bubbles, b, contact_map);

            // Unflip if original orientation was better
            if (flipped_score < score) {
                bubbles.flip(b);
            }
        }

        total_score = compute_total_consistency_score(bubbles, contact_map);

        phase_mutex.lock();
        if (total_score > best_score) {
            best_score = total_score;
            bubbles.get_phases(best_phases);
        }
        else {
            bubbles.set_phases(best_phases);
        }

        cerr << m << ' ' << best_score << ' ' << total_score << std::flush << '\n';
        phase_mutex.unlock();

        m = job_index.fetch_add(1);
    }
}


void phase_contacts(
        const contact_map_t& contact_map,
        const IncrementalIdMap<string>& id_map,
        Bubbles& bubbles,
        size_t n_threads
        ){

    vector<bool> best_phases(bubbles.size(), false);

    atomic<int64_t> best_score = std::numeric_limits<int64_t>::min();
    atomic<size_t> job_index = 0;

    size_t m_iterations = 10000;

    // Thread-related variables
    vector<thread> threads;
    mutex phase_mutex;

    // Launch threads
    for (uint64_t i=0; i<n_threads; i++){
        try {
            threads.emplace_back(thread(
                    random_phase_search,
                    cref(contact_map),
                    cref(id_map),
                    bubbles,
                    ref(best_phases),
                    ref(best_score),
                    ref(job_index),
                    ref(phase_mutex),
                    m_iterations
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

    bubbles.set_phases(best_phases);
}


void write_contact_map(
        path output_path,
        contact_map_t& contact_map,
        IncrementalIdMap<string>& id_map){
    ofstream output_file(output_path);

    if (not output_file.is_open() or not output_file.good()){
        throw std::runtime_error("ERROR: could not write to file: " + output_path.string());
    }

    for (auto& [id,map2]: contact_map){
        for (auto& [id2,count]: map2){
            output_file << id_map.get_name(id) << ',' << id_map.get_name(id2) << ',' << count << '\n';
        }
    }

}


void phase_hic(path sam_path, string required_prefix, int8_t min_mapq, size_t n_threads){
    // Id-to-name bimap for reference contigs
    IncrementalIdMap<string> id_map(true);

    // Mappings grouped by their read name (to simplify paired-end parsing)
    unpaired_mappings_t mappings;

    // Datastructures to represent linkages from hiC
    contact_map_t contact_map;
    vector <vector <int32_t> > adjacency;

    if (sam_path.extension() == ".sam") {
        parse_unpaired_sam_file(sam_path, mappings, id_map, required_prefix, min_mapq);
    }
    else if (sam_path.extension() == ".bam"){
        parse_unpaired_bam_file(sam_path, mappings, id_map, required_prefix, min_mapq);
    }
    else{
        throw runtime_error("ERROR: unrecognized extension for SAM/BAM input file: " + sam_path.extension().string());
    }

    remove_singleton_reads(mappings);

    // Build the contact map by iterating the pairs and creating edges in an all-by-all fashion between pairs
    generate_contact_map_from_mappings(mappings, id_map, contact_map);

    // To keep track of pairs of segments which exist in diploid bubbles
    Bubbles bubbles;

    // Initialize bubble objects using the shasta convention for bubbles
    generate_bubbles_from_shasta_names(bubbles, id_map);

    // Simplify contacts into bubble contacts (sometimes contacts are on one phase only)
    bubbles.generate_bubble_adjacency_from_contact_map(contact_map);

    cerr << "Phasing " << bubbles.size() << " bubbles" << '\n';

    print_mappings(mappings);

    generate_adjacency_matrix(bubbles, contact_map, adjacency);

    phase_contacts(contact_map, id_map, bubbles, n_threads);

    int64_t score = compute_total_consistency_score(bubbles, contact_map);

    path contacts_output_path = sam_path;
    path bandage_output_path = sam_path;

    string suffix1 = "p" + required_prefix;
    string suffix2 = "m" + to_string(int(min_mapq));
    string suffix3 = "s" + to_string(int(score));

    string contacts_suffix = suffix1 + "_" + suffix2 + "_" + suffix3 + "_contacts.csv";
    string bandage_suffix = suffix1 + "_" + suffix2 + "_" + suffix3 + "_bandage.csv";

    contacts_output_path.replace_extension(contacts_suffix);
    bandage_output_path.replace_extension(bandage_suffix);

    write_contact_map(contacts_output_path, contact_map, id_map);
    bubbles.write_bandage_csv(bandage_output_path, id_map);
}


int main (int argc, char* argv[]){
    path sam_path;
    string required_prefix;
    int8_t min_mapq = 0;
    size_t n_threads = 1;

    CLI::App app{"App description"};

    app.add_option(
            "-i,--input",
            sam_path,
            "Path to SAM or BAM containing filtered, paired HiC reads")
            ->required();

    app.add_option(
            "-p,--prefix",
            required_prefix,
            "Prefix required in ref name for mapping to be counted");

    app.add_option(
            "-m,--min_mapq",
            min_mapq,
            "Minimum required mapq value for mapping to be counted");

    app.add_option(
            "-t,--threads",
            n_threads,
            "Maximum number of threads to use");

    CLI11_PARSE(app, argc, argv);

    phase_hic(sam_path, required_prefix, min_mapq, n_threads);

    return 0;
}
