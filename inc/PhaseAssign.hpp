#ifndef GFASE_PHASEASSIGN_HPP
#define GFASE_PHASEASSIGN_HPP

#include "IncrementalIdMap.hpp"
#include "Filesystem.hpp"
#include "sparsepp/spp.h"
#include "Color.hpp"
#include "misc.hpp"
#include "Sam.hpp"

using gfase::rgb_to_hex;


namespace gfase{

class CigarSummary{
public:
    string primary_ref;
    uint32_t n_matches;
    uint32_t n_mismatches;
    uint32_t n_deletes;
    uint32_t n_inserts;
    uint32_t ref_length;
    uint32_t query_length;

    CigarSummary();
    void update(char cigar_operation, uint32_t length, uint32_t max_indel);
    double get_identity();
};


void assign_phases(
        path output_dir,
        path pat_ref_path,
        path mat_ref_path,
        path query_vs_pat_bam,
        path query_vs_mat_bam,
        path query_path,
        string required_prefix,
        size_t n_threads,
        unordered_map <string, array<CigarSummary,2> >& phased_cigar_summaries,
        array <set <string>, 2>& phased_contigs,
        map<string,size_t>& query_lengths,
        bool extract_fasta=false
);


void evaluate_phasing(
        path output_dir,
        path contact_phase_csv,
        path pat_ref_path,
        path mat_ref_path,
        path pat_bam_path,
        path mat_bam_path,
        path query_path,
        string required_prefix,
        size_t n_threads,
        bool extract_fasta
);

}


#endif //GFASE_PHASEASSIGN_HPP
