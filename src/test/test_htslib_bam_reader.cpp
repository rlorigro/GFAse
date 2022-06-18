#include "Filesystem.hpp"
#include "Bam.hpp"

using ghc::filesystem::path;
using gfase::SamElement;
using gfase::Bam;

#include <stdexcept>
#include <iostream>
#include <string>

using std::runtime_error;
using std::cerr;
using std::string;


int main(){
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();

    // Get test VCF path
    path relative_bam_path = "data/test.bam";
    path bam_path = project_directory / relative_bam_path;

    Bam bam_reader(bam_path);

    bam_reader.for_alignment_in_bam(true, [&](SamElement& e){
        cerr << e.ref_name << ' ' << e.query_name << ' ' << int(e.mapq) << ' ' << e.flag << '\n';

        e.for_each_cigar([&](auto type, auto length){
            cerr << type << length << ',';
        });
        cerr << '\n';
    });

    return 0;
}
