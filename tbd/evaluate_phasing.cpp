#include "bamtools/api/BamReader.h"
#include "IncrementalIdMap.hpp"
#include "Filesystem.hpp"
#include "sparsepp/spp.h"
#include "Color.hpp"
#include "CLI11.hpp"
#include "misc.hpp"
#include "Sam.hpp"

using gfase::for_element_in_sam_file;
using gfase::IncrementalIdMap;
using gfase::SamElement;
using gfase::rgb_to_hex;

using BamTools::BamAlignment;
using BamTools::BamReader;

using ghc::filesystem::create_directories;
using ghc::filesystem::exists;
using ghc::filesystem::path;
using spp::sparse_hash_set;
using CLI::App;

#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <stdexcept>
#include <fstream>
#include <utility>
#include <limits>
#include <vector>
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
using std::string;
using std::vector;
using std::array;
using std::pair;
using std::stoi;
using std::cerr;
using std::cref;
using std::ref;
using std::set;



void assign_phase(const string& ref_name, const string& read_name, array <set <string>, 2>& phases){
    bool a = ref_name.find("pat") != std::string::npos;
    bool b = ref_name.find("mat") != std::string::npos;

    if (a and !b){
//        cerr << ref_name << 0 << '\n';
        phases[0].emplace(read_name);
    }
    else if (!a and b){
//        cerr << ref_name << 1 << '\n';
        phases[1].emplace(read_name);
    }
    else if (a and b){
        throw std::runtime_error("ERROR: found mat and pat label in single ref contig: " + ref_name);
    }

}


void get_reference_lengths_from_bam_file(path bam_path, unordered_map<string,size_t>& contig_lengths){
    BamReader reader;

    if (!reader.Open(bam_path) ) {
        throw std::runtime_error("ERROR: could not read BAM file: " + bam_path.string());
    }

    auto& ref_data = reader.GetReferenceData();

    for (auto& item: ref_data){
        contig_lengths.emplace(item.RefName, item.RefLength);
    }
}


void get_query_lengths_from_bam_file(path bam_path, unordered_map<string,size_t>& contig_lengths){
    BamReader reader;

    if (!reader.Open(bam_path) ) {
        throw std::runtime_error("ERROR: could not read BAM file: " + bam_path.string());
    }

    auto& ref_data = reader.GetReferenceData();

    for (auto& item: ref_data){
        contig_lengths.emplace(item.RefName, item.RefLength);
    }
}


void parse_unpaired_bam_file(
        path bam_path,
        array <set <string>, 2>& true_phases,
        array <set <string>, 2>& inferred_phases,
        unordered_map<string,size_t>& query_lengths
        ){

    BamReader reader;

    if (!reader.Open(bam_path) ) {
        throw std::runtime_error("ERROR: could not read BAM file: " + bam_path.string());
    }

    auto reference_data = reader.GetReferenceData();

    BamAlignment e;

    size_t l = 0;
    while (reader.GetNextAlignment(e) ) {
        if (e.IsPrimaryAlignment() and (not e.IsSupplementary()) and e.IsMapped()){
            auto& ref_name = reference_data.at(e.RefID).RefName;

            // Only add entries that are in both sets
            bool found = ((inferred_phases[0].count(e.Name) > 0) or (inferred_phases[1].count(e.Name) > 0));

            if (found) {
//                cerr << l << ' ' << e.Name << ' ' << ref_name << ' ' << e.MapQuality << ' ' << e.IsPrimaryAlignment() << ' ' << e.IsMapped() <<'\n';
                assign_phase(ref_name, e.Name, true_phases);
                query_lengths.emplace(e.Name, e.Length);
            }
        }
        l++;
    }
}


void for_entry_in_csv(path csv_path, const function<void(string& name, bool phase)>& f){
    ifstream file(csv_path);

    if (not file.is_open() or not file.good()){
        throw runtime_error("ERROR: could not read input file: " + csv_path.string());
    }

    char c;

    size_t n_delimiters = 0;
    int64_t n_lines = 0;

    string ref_name;
    string phase_token;

    file.ignore(numeric_limits<streamsize>::max(), '\n');

    while (file.get(c)){
        if (c == '\n'){
            n_delimiters = 0;

//            cerr << n_lines << ' ' << ref_name << ' ' << phase_token << '\n';

            f(ref_name, bool(stoi(phase_token)));

            ref_name.clear();
            phase_token.clear();

            n_lines++;
        }
        else if (c == ','){
            n_delimiters++;
        }
        else{
            if (n_delimiters == 0){
                ref_name += c;
            }
            else if (n_delimiters == 1){
                phase_token += c;
            }
        }
    }
}


void parse_phase_csv(path csv_path, array <set <string>, 2>& phases){
    for_entry_in_csv(csv_path, [&](string& name, bool phase){
        phases[phase].emplace(name);
    });
}


/*
    #####  Color Palette by Paletton.com
    #####  Palette URL: http://paletton.com/#uid=7000u0kr6rGgDzJlztwt9l+yoh7

    *** Primary color:

       shade 0 = #DD2222 = rgb(221, 34, 34) = rgba(221, 34, 34,1) = rgb0(0.867,0.133,0.133)

    *** Secondary color (1):

       shade 0 = #DD7622 = rgb(221,118, 34) = rgba(221,118, 34,1) = rgb0(0.867,0.463,0.133)

    *** Secondary color (2):

       shade 0 = #148484 = rgb( 20,132,132) = rgba( 20,132,132,1) = rgb0(0.078,0.518,0.518)

    *** Complement color:

       shade 0 = #1BB01B = rgb( 27,176, 27) = rgba( 27,176, 27,1) = rgb0(0.106,0.69,0.106)

    #####  Generated by Paletton.com (c) 2002-2014
 */
void write_bandage_csv(
        const vector<string>& intersection_00,
        const vector<string>& intersection_11,
        const vector<string>& intersection_01,
        const vector<string>& intersection_10,
        path output_path){

    ofstream file(output_path);

    file << "Name" << ',' << "Set" << ',' << "Color" << '\n';
    for (auto& item: intersection_00){
        file << item << ',' << "00" << ',' << '#' << rgb_to_hex(0.867,0.133,0.133) << '\n';     // Red
    }
    for (auto& item: intersection_11){
        file << item << ',' << "11" << ',' << '#' << rgb_to_hex(0.867,0.463,0.133) << '\n';     // Orange
    }
    for (auto& item: intersection_01){
        file << item << ',' << "01" << ',' << '#' << rgb_to_hex(0.078,0.518,0.518) << '\n';     // Blue
    }
    for (auto& item: intersection_10){
        file << item << ',' << "10" << ',' << '#' << rgb_to_hex(0.106,0.69,0.106) << '\n';      // Green
    }
}


size_t get_total_length_of_phase_set(const unordered_map<string,size_t>& contig_lengths, const vector<string>& s){
    size_t l = 0;

    for (auto& item: s){
        l += contig_lengths.at(item);
    }

    return l;
}


void evaluate_phasing(path sam_path, path csv_path, path output_dir){
    if (exists(output_dir)){
        throw runtime_error("ERROR: output directory must not exist already, remove it or choose another directory: " + output_dir.string());
    }

    create_directories(output_dir);

    array <set <string>, 2> true_phases;
    array <set <string>, 2> inferred_phases;
    unordered_map<string,size_t> query_lengths; // TODO: FIX THIS, CURRENTLY DOES NOT FIND TRUE QUERY LENGTH

    parse_phase_csv(csv_path, inferred_phases);

    if (sam_path.extension() == ".bam"){
        parse_unpaired_bam_file(sam_path, true_phases, inferred_phases, query_lengths);
    }
    else{
        throw runtime_error("ERROR: unrecognized extension for BAM input file: " + sam_path.extension().string());
    }

//    cerr << "Inferred Phase 0" << '\n';
//    for (auto& item: inferred_phases[0]){
//        cerr << '\t' << item << '\n';
//    }
//    cerr << "Inferred Phase 1" << '\n';
//    for (auto& item: inferred_phases[1]){
//        cerr << '\t' << item << '\n';
//    }
//
//    cerr << "True Phase 0" << '\n';
//    for (auto& item: true_phases[0]){
//        cerr << '\t' << item << '\n';
//    }
//    cerr << "True Phase 1" << '\n';
//    for (auto& item: true_phases[1]){
//        cerr << '\t' << item << '\n';
//    }

    for (auto& item: query_lengths){
        cerr << item.first << ' ' << item.second << '\n';
    }

    vector<string> intersection_00;
    vector<string> intersection_11;
    vector<string> intersection_01;
    vector<string> intersection_10;

    set_intersection(true_phases[0].begin(), true_phases[0].end(),
                     inferred_phases[0].begin(), inferred_phases[0].end(),
                     std::back_inserter(intersection_00));

    set_intersection(true_phases[1].begin(), true_phases[1].end(),
                     inferred_phases[1].begin(), inferred_phases[1].end(),
                     std::back_inserter(intersection_11));

    set_intersection(true_phases[0].begin(), true_phases[0].end(),
                     inferred_phases[1].begin(), inferred_phases[1].end(),
                     std::back_inserter(intersection_01));

    set_intersection(true_phases[1].begin(), true_phases[1].end(),
                     inferred_phases[0].begin(), inferred_phases[0].end(),
                     std::back_inserter(intersection_10));

    cerr << "Intersection" << '\n';
    cerr << '\t' << "00" << ' ' << intersection_00.size() << ' ' << get_total_length_of_phase_set(query_lengths, intersection_00) << '\n';
    cerr << '\t' << "11" << ' ' << intersection_11.size() << ' ' << get_total_length_of_phase_set(query_lengths, intersection_11) << '\n';
    cerr << '\t' << "01" << ' ' << intersection_01.size() << ' ' << get_total_length_of_phase_set(query_lengths, intersection_01) << '\n';
    cerr << '\t' << "10" << ' ' << intersection_10.size() << ' ' << get_total_length_of_phase_set(query_lengths, intersection_10) << '\n';

    write_bandage_csv(intersection_00, intersection_11, intersection_01, intersection_10, output_dir / "phase_intersections.csv");
}


int main (int argc, char* argv[]){
    path output_dir;
    path bam_path;
    path phase_csv;

    CLI::App app{"App description"};

    app.add_option(
            "-o,--output_dir",
            output_dir,
            "Path to directory which will be created and contain the results of the evaluation. Directory must not exist yet.")
            ->required();

    app.add_option(
            "-r,--reference",
            bam_path,
            "Path to BAM containing alignments of assembly contigs to a phased reference. Reference contigs must contain 'mat' or 'pat' in their name.")
            ->required();

    app.add_option(
            "-c,--csv",
            phase_csv,
            "Path to CSV file containing fields Name,Phase at columns 0,1")
            ->required();


    CLI11_PARSE(app, argc, argv);

    evaluate_phasing(bam_path, phase_csv, output_dir);

    return 0;
}