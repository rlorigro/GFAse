#include "Filesystem.hpp"
#include "CLI11.hpp"
#include "Sam.hpp"

using gfase::for_element_in_sam_file;
using ghc::filesystem::path;
using gfase::SamElement;
using CLI::App;

#include <stdexcept>
#include <ostream>
#include <string>
#include <map>

using std::runtime_error;
using std::ifstream;
using std::ofstream;
using std::string;
using std::map;


void get_mapq_distribution(path sam_path){
    map<int8_t,size_t> distribution;

    for_element_in_sam_file(sam_path, [&](SamElement& e){
        if (not e.is_not_primary()) {
            distribution[e.mapq]++;
        }
    });

    path output_path = sam_path;
    output_path.replace_extension("mapq_distribution.csv");

    ofstream file(output_path);

    if ((not file.is_open()) or (not file.good())){
        throw runtime_error("ERROR: file could not be written: " + output_path.string());
    }

    for (auto& [key, frequency]: distribution){
        file << int(key) << ',' << frequency << '\n';
    }
}


int main (int argc, char* argv[]){
    path sam_path;

    CLI::App app{"App description"};

    app.add_option(
            "-i,--input",
            sam_path,
            "Path to SAM containing filtered, paired HiC reads")
            ->required();

    CLI11_PARSE(app, argc, argv);

    get_mapq_distribution(sam_path);

    return 0;
}

