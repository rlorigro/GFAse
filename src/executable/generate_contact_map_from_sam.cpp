#include "IncrementalIdMap.hpp"
#include "SAMElement.hpp"
#include "Filesystem.hpp"
#include "CLI11.hpp"

using ghc::filesystem::path;
using gfase::IncrementalIdMap;
using gfase::sam_comparator;
using gfase::SAMElement;
using CLI::App;

#include <unordered_map>
#include <unordered_set>
#include <stdexcept>
#include <fstream>
#include <utility>
#include <random>
#include <limits>
#include <bitset>
#include <vector>
#include <array>
#include <set>

using std::unordered_set;
using std::unordered_map;
using std::runtime_error;
using std::numeric_limits;
using std::streamsize;
using std::to_string;
using std::ifstream;
using std::ofstream;
using std::string;
using std::vector;
using std::bitset;
using std::array;
using std::pair;
using std::stoi;
using std::cerr;
using std::set;


using mappings_t = unordered_map <string, array <set <SAMElement, decltype(sam_comparator)*>, 2> >;


void print_mappings(const mappings_t& mappings){
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


class Bubble{
public:
    array<int32_t,2> ids;
    bool phase;

    Bubble();
    Bubble(int32_t id1, int32_t id2, bool phase);
    void flip();
    int32_t first();
    int32_t second();
    int32_t get(bool i);
    int32_t is_first(int32_t id);
    int32_t is_second(int32_t id);
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


int32_t Bubble::first(){
    return ids[0 + phase];
}


int32_t Bubble::second(){
    return ids[1 - phase];
}

///
/// s p return
/// 0 1 1
/// 1 0 1
/// 0 0 0
/// 1 1 0
int32_t Bubble::get(bool side){
    return ids[side != phase];
}


int32_t Bubble::is_first(int32_t id){
    if (get(id) == 0){
        return true;
    }
    else{
        return false;
    }
}


int32_t Bubble::is_second(int32_t id){
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
    unordered_map<int32_t,int32_t> id_to_bubble;

    Bubbles();
    void write_bandage_csv(path output_path, IncrementalIdMap<string>& id_map);
    void at(size_t i, Bubble& b);
    void emplace(int32_t id1, int32_t id2, bool phase);
    int32_t find(int32_t id);
    void flip(size_t b);
    size_t size() const;
};


Bubbles::Bubbles():
    bubbles(),
    id_to_bubble()
{}


void Bubbles::write_bandage_csv(path output_path, IncrementalIdMap<string>& id_map){
    ofstream file(output_path);

    if (not file.is_open() or not file.good()){
        throw std::runtime_error("ERROR: could not write to file: " + output_path.string());
    }

    file << "Name" << ',' << "Color" << '\n';

    for (auto& b: bubbles){
        auto id0 = b.get(0);
        auto id1 = b.get(1);

        file << id_map.get_name(id0) << ',' << "Dark Orange" << '\n';
        file << id_map.get_name(id1) << ',' << "Green Yellow" << '\n';
    }


}


void Bubbles::emplace(int32_t id1, int32_t id2, bool phase){
    id_to_bubble[id1] = int32_t(bubbles.size());
    id_to_bubble[id2] = int32_t(bubbles.size());
    bubbles.emplace_back(id1, id2, phase);
}


size_t Bubbles::size() const{
    return bubbles.size();
}


void Bubbles::at(size_t i, Bubble& b){
    b = bubbles.at(i);
}


int32_t Bubbles::find(int32_t id){
    return id_to_bubble.at(id);
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

        cerr << name << ' ' << side << ' ' << other_name << '\n';

        visited.emplace(id);
        visited.emplace(other_id);
    }
}


void parse_sam_file(
        path sam_path,
        mappings_t& mappings,
        IncrementalIdMap<string>& id_map,
        string required_prefix,
        int8_t min_mapq){

    ifstream file(sam_path);

    if (not file.is_open() or not file.good()){
        throw runtime_error("ERROR: could not read input file: " + sam_path.string());
    }

    char c;
    size_t n_delimiters = 0;
    int32_t n_lines = 0;

    string read_name;
    string flag_token;
    string ref_name;
    string mapq_token;

    char header_delimiter = '@';

    // If this is a header line, skip to the next line and increment n_lines until we are no longer on a header
    while (file.peek() == header_delimiter){
        cerr << string(1, file.peek()) << '\n';
        file.ignore(numeric_limits<streamsize>::max(), '\n');
        n_lines++;
    }

    while (file.get(c)){
        if (c == '\n'){
            n_delimiters = 0;

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
                auto mapq = int8_t(stoi(mapq_token));

                if (mapq >= min_mapq) {
                    SAMElement e(read_name, ref_name, n_lines, stoi(flag_token), mapq);
                    mappings[read_name][e.is_second_mate()].emplace(e);
                }
            }

//            cerr << n_lines << " " << read_name << " " << ref_name << '\n';

            read_name.clear();
            flag_token.clear();
            ref_name.clear();
            mapq_token.clear();

            n_lines++;
        }
        else if (c == '\t'){
            n_delimiters++;
        }
        else{
            if (n_delimiters == 0){
                read_name += c;
            }
            else if (n_delimiters == 1){
                flag_token += c;
            }
            else if (n_delimiters == 2){
                ref_name += c;
            }
            else if (n_delimiters == 4){
                mapq_token += c;
            }
        }
    }
}


void remove_unpaired_reads(mappings_t& mappings){
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


void generate_contact_map_from_mappings(
        mappings_t& mappings,
        IncrementalIdMap<string>& id_map,
        unordered_map <int32_t, unordered_map<int32_t, int32_t> >& contact_map
){
    for (auto& [name,mates]: mappings) {
        size_t i = 0;

        for (auto& e: mates[0]){
            size_t j = 0;

            for (auto& e2: mates[1]){
                if (j >= i){
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


void generate_adjacency_matrix(
        const Bubbles& bubbles,
        const unordered_map <int32_t,unordered_map<int32_t, int32_t> >& contact_map,
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


int64_t compute_consistency_score(
        Bubbles& bubbles,
        size_t bubble_index,
        const unordered_map <int32_t, unordered_map<int32_t, int32_t> >& contact_map
        ){
    int64_t score = 0;

    Bubble bubble;
    bubbles.at(int32_t(bubble_index), bubble);

    auto id0 = bubble.get(0);
    auto id1 = bubble.get(1);

    unordered_set<size_t> other_bubbles;

    auto c0 = contact_map.find(id0);
    auto c1 = contact_map.find(id1);

    // Iterate all contacts with id0
    if (c0 != contact_map.end()) {
        for (auto&[other_id, count]: c0->second) {
            Bubble other_bubble;
            auto b_other = bubbles.find(other_id);

            other_bubbles.emplace(b_other);
        }
    }

    // Iterate all contacts with id1
    if (c1 != contact_map.end()) {
        for (auto&[other_id, count]: c1->second) {
            Bubble other_bubble;
            auto b_other = bubbles.find(other_id);

            other_bubbles.emplace(b_other);
        }
    }

    for (auto b_other: other_bubbles){
        Bubble other_bubble;
        bubbles.at(b_other, other_bubble);

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
    }

    return score;
}


void phase_contacts(
        const unordered_map <int32_t, unordered_map<int32_t, int32_t> >& contact_map,
        const IncrementalIdMap<string>& id_map,
        Bubbles& bubbles){

    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<int> uniform_distribution(0,int(bubbles.size()-1));

    int64_t prev_total_score = std::numeric_limits<int64_t>::min();
    int64_t total_score = 0;

    for (size_t m=0; m<10; m++) {
        for (size_t i=0; i<20; i++) {
//            cerr << "--- " << i << " ---" << '\n';
            for (size_t j=0; j<bubbles.size(); j++) {
//                auto b = j;
                auto b = uniform_distribution(rng);

                auto score = compute_consistency_score(bubbles, b, contact_map);
                bubbles.flip(b);
                auto flipped_score = compute_consistency_score(bubbles, b, contact_map);

                // Unflip if original orientation was better
                if (flipped_score < score) {
                    bubbles.flip(b);
                }

//                cerr << b << ' ' << score << ' ' << flipped_score << '\n';
            }
        }

        total_score = 0;
        for (size_t b=0; b<bubbles.size(); b++){
            total_score += compute_consistency_score(bubbles, b, contact_map);
        }
        prev_total_score = total_score;
        cerr << m << ' ' << total_score << '\n';

        // Randomly perturb
        for (size_t i=0; i<((bubbles.size()/100) + 1); i++) {
            bubbles.flip(uniform_distribution(rng));
        }
    }
}


void write_contact_map(
        path output_path,
        unordered_map <int32_t, unordered_map<int32_t, int32_t> >& contact_map,
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


void phase_hic(path sam_path, string required_prefix, int8_t min_mapq){
    // Id-to-name bimap for reference contigs
    IncrementalIdMap<string> id_map(true);

    // Mappings grouped by their read name (to simplify paired-end parsing)
    unordered_map <string, array <set <SAMElement, decltype(sam_comparator)*>, 2> > mappings;

    // Datastructures to represent linkages from hiC
    unordered_map <int32_t, unordered_map<int32_t, int32_t> > contact_map;
    vector <vector <int32_t> > adjacency;

    // To keep track of pairs of segments which exist in diploid bubbles
    Bubbles bubbles;

    parse_sam_file(sam_path, mappings, id_map, required_prefix, min_mapq);

    // Filter the mappings for each read pair to verify that they actually have a first and second pair
    remove_unpaired_reads(mappings);

    // Build the contact map by iterating the pairs and creating edges in an all-by-all fashion between pairs
    generate_contact_map_from_mappings(mappings, id_map, contact_map);

    // Initialize bubble objects using the shasta convention for bubbles
    generate_bubbles_from_shasta_names(bubbles, id_map);

//    print_mappings(mappings);

    generate_adjacency_matrix(bubbles, contact_map, adjacency);

    phase_contacts(contact_map, id_map, bubbles);

    path contacts_output_path = sam_path;
    path bandage_output_path = sam_path;

    string suffix1 = "p" + required_prefix;
    string suffix2 = "m" + to_string(int(min_mapq));

    string contacts_suffix = suffix1 + "_" + suffix2 + "_contacts.csv";
    string bandage_suffix = suffix1 + "_" + suffix2 + "_bandage.csv";

    contacts_output_path.replace_extension(contacts_suffix);
    bandage_output_path.replace_extension(bandage_suffix);

    write_contact_map(contacts_output_path, contact_map, id_map);
    bubbles.write_bandage_csv(bandage_output_path, id_map);
}


int main (int argc, char* argv[]){
    path sam_path;
    string required_prefix;
    int8_t min_mapq = 0;

    CLI::App app{"App description"};

    app.add_option(
            "-i,--input",
            sam_path,
            "Path to SAM containing filtered, paired HiC reads")
            ->required();

    app.add_option(
            "-p,--prefix",
            required_prefix,
            "Prefix required in ref name for mapping to be counted");

    app.add_option(
            "-m,--min_mapq",
            min_mapq,
            "Minimum required mapq value for mapping to be counted");

    CLI11_PARSE(app, argc, argv);

    phase_hic(sam_path, required_prefix, min_mapq);

    return 0;
}
