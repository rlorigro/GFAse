#ifndef GFASE_BAM_HPP
#define GFASE_BAM_HPP

#include "htslib/include/htslib/hts.h"
#include "htslib/include/htslib/sam.h"
#include "MultiContactGraph.hpp"
#include "Filesystem.hpp"
#include "Sam.hpp"

using ghc::filesystem::path;

#include <functional>
#include <string>

using std::function;
using std::string;


namespace gfase {

class Bam {
    path bam_path;

    samFile* bam_file;
    bam_hdr_t* bam_header;
    hts_itr_t* bam_iterator;
    bam1_t* alignment;

public:
    Bam(path bam_path);
    ~Bam();
    void for_alignment_in_bam(const function<void(const string& ref_name, const string& query_name, uint8_t map_quality, uint16_t flag)>& f);
    void for_alignment_in_bam(bool get_cigar, const function<void(SamElement& alignment)>& f);
    void for_alignment_in_bam(const function<void(FullAlignmentBlock& a)>& f);
    void for_ref_in_header(const function<void(const string& ref_name, uint32_t length)>& f) const;
    static bool is_first_mate(uint16_t flag);
    static bool is_second_mate(uint16_t flag);
    static bool is_not_primary(uint16_t flag);
    static bool is_primary(uint16_t flag);
    static bool is_supplementary(uint16_t flag);
};


template <class T> void update_contact_map(
        const vector<T>& alignments,
        MultiContactGraph& contact_graph,
        IncrementalIdMap<string>& id_map);


void parse_unpaired_bam_file(
        path bam_path,
        MultiContactGraph& contact_graph,
        IncrementalIdMap<string>& id_map,
        string required_prefix,
        int8_t min_mapq);


}

#endif //GFASE_BAM_HPP
