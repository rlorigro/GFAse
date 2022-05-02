#include "OverlapMap.hpp"
#include "PafElement.hpp"

#include "boost/program_options.hpp"
#include "boost/bimap.hpp"

using boost::program_options::options_description;
using boost::program_options::variables_map;
using boost::program_options::value;
using boost::bimap;


#include "Filesystem.hpp"
#include <iostream>
#include <fstream>
#include <utility>
#include <string>
#include <vector>

using ghc::filesystem::create_directories;
using ghc::filesystem::path;
using std::runtime_error;
using std::ifstream;
using std::ofstream;
using std::to_string;
using std::string;
using std::vector;
using std::pair;
using std::cerr;
using std::cout;

typedef bimap<uint32_t,string> uint32_string_bimap;
typedef uint32_string_bimap::value_type bimap_pair;

using gfase::RegionalOverlapMap;


/// Parse a PAF file: https://github.com/lh3/miniasm/blob/master/PAF.md using the following data:
///
/// Col (1-based)       Data                        Type
/// -------------       ----                        ----
/// 1                   "Query sequence name"       string (assumed to be numeric for this project, using ONT reads)
/// 6                   "Target sequence name"      string
/// 8                   "Target start..."           int
/// 9                   "Target end..."             int
/// 12                  "Mapping quality"           int
///
void load_paf_as_overlap_map(
        path paf_path,
        uint32_string_bimap& id_vs_name,
        RegionalOverlapMap& overlap_map,
        uint32_t min_quality) {

    ifstream paf_file(paf_path);

    if (not paf_file.good()) {
        throw runtime_error("ERROR: could not open input file: " + paf_path.string());
    }

    string token;
    string region_name;
    string read_name;
    uint32_t start = 0;
    uint32_t stop = 0;
    uint32_t quality = 0;

    uint64_t n_delimiters = 0;
    uint64_t n_lines = 0;
    char c;

    while (paf_file.get(c)) {
        if (c == '\t') {
            if (n_delimiters == 0) {
                read_name = token;
            } else if (n_delimiters == 5) {
                region_name = token;
            } else if (n_delimiters == 7) {
                start = stoi(token);
            } else if (n_delimiters == 8) {
                stop = stoi(token);
            } else if (n_delimiters == 11) {
                quality = stoi(token);

                if (quality >= min_quality) {
                    auto result = id_vs_name.right.find(read_name);
                    uint32_t id;

                    // The bimap for id <-> name is not necessarily initialized for this read.
                    // If it does exist, then just fetch the ID
                    if (result != id_vs_name.right.end()) {
                        id = result->second;
                    }
                        // If it doesn't exist, then make a new ID by incrementing by 1 (aka get the size)
                    else {
                        id = id_vs_name.size();
                        id_vs_name.insert(bimap_pair(id, read_name));
                    }

                    overlap_map.insert(region_name, start, stop, id);
                }
            }

            token.resize(0);
            n_delimiters++;
        } else if (c == '\n') {
            if (n_delimiters < 11) {
                throw runtime_error(
                        "ERROR: file provided does not contain sufficient tab delimiters to be PAF on line: " +
                        to_string(n_lines));
            }

            token.resize(0);
            n_delimiters = 0;
            n_lines++;
        } else {
            token += c;
        }
    }
}


void find_overlapping_alignments(path paf_path) {
    uint32_string_bimap id_vs_name;
    RegionalOverlapMap overlap_map;
    uint32_t min_quality = 10;

    load_paf_as_overlap_map(paf_path, id_vs_name, overlap_map, min_quality);

    for (auto& item: id_vs_name) {
        cout << item.get_left() << ',' << item.get_right() << '\n';
    }
    cout << '\n';
    overlap_map.print(cout);
}



int main(int argc, char* argv[]) {
    path paf_path;

    options_description options("Arguments:");

    options.add_options()
            ("paf_path",
             value<path>(&paf_path)->required(),
             "File path of PAF file containing alignments to some reference");

    variables_map vm;
    store(parse_command_line(argc, argv, options), vm);

    // If help was specified, or no arguments given, provide help
    if (vm.count("help") || argc == 1) {
        cerr << options << "\n";
        return 0;
    }
    notify(vm);

    find_overlapping_alignments(paf_path);

    return 0;
}
