#include "Filesystem.hpp"
#include "CLI11.hpp"
#include "Bam.hpp"
#include "Sam.hpp"

using gfase::for_element_in_sam_file;
using gfase::SamElement;
using gfase::Bam;

using ghc::filesystem::path;
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

    if (sam_path.extension() == ".sam") {
        for_element_in_sam_file(sam_path, [&](SamElement& e) {
            if (not e.is_not_primary()) {
                distribution[e.mapq]++;
            }
        });
    }
    else if (sam_path.extension() == ".bam"){
        Bam reader(sam_path);

        size_t l = 0;
        reader.for_alignment_in_bam(false, [&](SamElement& a) {
            distribution[int8_t(a.mapq)]++;
            l++;
        });
    }
    else{
        throw runtime_error("ERROR: file format not bam or sam: " + sam_path.string());
    }

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

