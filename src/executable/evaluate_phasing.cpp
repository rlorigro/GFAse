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



void assign_phase(const string& ref_name, const string& read_name, array <sparse_hash_set <string>, 2>& phases){
    bool a = ref_name.find("pat") != std::string::npos;
    bool b = ref_name.find("mat") != std::string::npos;

    if (a and !b){
        cerr << ref_name << 0 << '\n';
        phases[0].emplace(read_name);
    }
    else if (!a and b){
        cerr << ref_name << 1 << '\n';
        phases[1].emplace(read_name);
    }
    else if (a and b){
        throw std::runtime_error("ERROR: found mat and pat label in single ref contig: " + ref_name);
    }

}


void parse_unpaired_sam_file(
        path sam_path,
        array <sparse_hash_set <string>, 2>& true_phases,
        array <sparse_hash_set <string>, 2>& inferred_phases
        ){

    for_element_in_sam_file(sam_path, [&](SamElement& e){
        if (not e.is_not_primary() and not e.is_supplementary()){
            // Only add entries that are in both sets
            bool found = ((inferred_phases[0].count(e.ref_name) > 0) or (inferred_phases[1].count(e.ref_name) > 0));

            if (found) {
                assign_phase(e.ref_name, e.read_name, true_phases);
            }
        }
    });
}


void parse_unpaired_bam_file(
        path bam_path,
        array <sparse_hash_set <string>, 2>& true_phases,
        array <sparse_hash_set <string>, 2>& inferred_phases
        ){

    BamReader reader;

    if (!reader.Open(bam_path) ) {
        throw std::runtime_error("ERROR: could not read BAM file: " + bam_path.string());
    }

    auto reference_data = reader.GetReferenceData();

    BamAlignment e;

    while (reader.GetNextAlignment(e) ) {
        auto& ref_name = reference_data[e.RefID].RefName;

        if (e.IsPrimaryAlignment() and not e.IsSupplementary()){
            // Only add entries that are in both sets
            bool found = ((inferred_phases[0].count(ref_name) > 0) or (inferred_phases[1].count(ref_name) > 0));

            if (found) {
                assign_phase(ref_name, e.Name, true_phases);
            }
        }
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

            cerr << n_lines << ' ' << ref_name << ' ' << phase_token << '\n';

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


void parse_phase_csv(path csv_path, array <sparse_hash_set <string>, 2>& phases){
    for_entry_in_csv(csv_path, [&](string& name, bool phase){
        phases[phase].emplace(name);
    });
}



void evaluate_phasing(path sam_path, path csv_path){
    array <sparse_hash_set <string>, 2> true_phases;
    array <sparse_hash_set <string>, 2> inferred_phases;

    parse_phase_csv(csv_path, inferred_phases);

    if (sam_path.extension() == ".sam") {
        parse_unpaired_sam_file(sam_path, true_phases, inferred_phases);
    }
    else if (sam_path.extension() == ".bam"){
        parse_unpaired_bam_file(sam_path, true_phases, inferred_phases);
    }
    else{
        throw runtime_error("ERROR: unrecognized extension for SAM/BAM input file: " + sam_path.extension().string());
    }

    vector<string> intersection_00;
    vector<string> intersection_11;
    vector<string> intersection_01;
    vector<string> intersection_10;

    std::set_intersection(true_phases[0].begin(), inferred_phases[0].end(),
                          true_phases[0].begin(), inferred_phases[0].end(),
                          std::back_inserter(intersection_00));


    std::set_intersection(true_phases[1].begin(), inferred_phases[1].end(),
                          true_phases[1].begin(), inferred_phases[1].end(),
                          std::back_inserter(intersection_11));


    std::set_intersection(true_phases[0].begin(), inferred_phases[1].end(),
                          true_phases[0].begin(), inferred_phases[1].end(),
                          std::back_inserter(intersection_01));


    std::set_intersection(true_phases[1].begin(), inferred_phases[0].end(),
                          true_phases[1].begin(), inferred_phases[0].end(),
                          std::back_inserter(intersection_10));


    cerr << "00" << intersection_00.size() << '\n';
    cerr << "11" << intersection_11.size() << '\n';
    cerr << "01" << intersection_01.size() << '\n';
    cerr << "10" << intersection_10.size() << '\n';

}


int main (int argc, char* argv[]){
    path sam_path;
    path phase_csv;

    CLI::App app{"App description"};

    app.add_option(
            "-r,--reference",
            sam_path,
            "Path to SAM or BAM containing alignments of assembly contigs to a phased reference. Reference contigs must contain 'mat' or 'pat' in their name.")
            ->required();

    app.add_option(
            "-c,--csv",
            phase_csv,
            "Path to CSV file containing fields Name,Phase at columns 0,1")
            ->required();


    CLI11_PARSE(app, argc, argv);

    evaluate_phasing(sam_path, phase_csv);

    return 0;
}