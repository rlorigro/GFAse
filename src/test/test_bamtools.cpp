#include "bamtools/api/BamMultiReader.h"
#include "Filesystem.hpp"

using BamTools::BamReader;
using BamTools::BamAlignment;
using ghc::filesystem::path;

#include <iostream>

using std::cerr;



int main(){
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();

    // Get test VCF path
    path relative_bam_path = "data/test.bam";
    path absolute_bam_path = project_directory / relative_bam_path;

    cerr << 0 << '\n' << std::flush;
    BamReader reader;

    cerr << 1 << '\n' << std::flush;

    if (!reader.Open(absolute_bam_path) ) {
        cerr << 2 << '\n' << std::flush;
        cerr << "Could not open input BAM files." << '\n';
    }

    cerr << 3 << '\n' << std::flush;

    BamAlignment a;
    while ( reader.GetNextAlignmentCore(a) ) {
        if ( a.MapQuality == 59 )
            cerr << '\t' << a.MapQuality << '\n' << std::flush;
    }

    cerr << 4 << '\n' << std::flush;

    reader.Close();

    cerr << 5 << '\n' << std::flush;

    return 0;
}
