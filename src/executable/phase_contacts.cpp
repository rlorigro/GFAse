#include "bamtools/api/BamReader.h"
#include "IncrementalIdMap.hpp"
#include "BubbleGraph.hpp"
#include "Filesystem.hpp"
#include "sparsepp/spp.h"
#include "CLI11.hpp"
#include "misc.hpp"
#include "Sam.hpp"

using gfase::for_element_in_sam_file;
using gfase::IncrementalIdMap;
using gfase::paired_mappings_t;
using gfase::BubbleGraph;
using gfase::unpaired_mappings_t;
using gfase::SamElement;
using gfase::contact_map_t;

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
    BubbleGraph bubble_graph(id_map, contact_map);

    cerr << "Phasing " << bubble_graph.size() << " bubbles" << '\n';

    print_mappings(mappings);

    generate_adjacency_matrix(bubble_graph, contact_map, adjacency);

    phase_contacts(contact_map, id_map, bubble_graph, n_threads);

    int64_t score = compute_total_consistency_score(bubble_graph, contact_map);

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
    bubble_graph.write_bandage_csv(bandage_output_path, id_map);
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
