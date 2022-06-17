#include "Filesystem.hpp"
#include "htslib/include/htslib/hts.h"
#include "htslib/include/htslib/sam.h"

using ghc::filesystem::path;

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

    samFile* bam_file;
    bam_hdr_t* bam_header;
    hts_idx_t* bam_index;
    hts_itr_t* bam_iterator;
    bam1_t* alignment;

    bam_file = nullptr;
    bam_index = nullptr;
    bam_iterator = nullptr;
    alignment = bam_init1();

    if ((bam_file = hts_open(bam_path.string().c_str(), "r")) == 0) {
        throw runtime_error("ERROR: Cannot open bam file: " + bam_path.string());
    }

    // bam index
    if ((bam_index = sam_index_load(bam_file, bam_path.string().c_str())) == 0) {
        throw runtime_error("ERROR: Cannot open index for bam file: " + bam_path.string() + "\n");
    }

    // bam header
    if ((bam_header = sam_hdr_read(bam_file)) == 0){
        throw runtime_error("ERROR: Cannot open header for bam file: " + bam_path.string() + "\n");
    }

    while (sam_read1(bam_file, bam_header, alignment) >= 0){
        string read_name = bam_get_qname(alignment);
        string ref_name = bam_header->target_name[alignment->core.tid];

        cerr << read_name << ' ' << ref_name << '\n';
    }

    hts_close(bam_file);
    bam_hdr_destroy(bam_header);
    bam_destroy1(alignment);
    hts_idx_destroy(bam_index);
    hts_itr_destroy(bam_iterator);

    return 0;
}