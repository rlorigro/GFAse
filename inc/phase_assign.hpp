#ifndef GFASE_PHASE_ASSIGN_HPP
#define GFASE_PHASE_ASSIGN_HPP

#include "bamtools/api/BamReader.h"
#include "IncrementalIdMap.hpp"
#include "Filesystem.hpp"
#include "sparsepp/spp.h"
#include "Color.hpp"
#include "misc.hpp"
#include "Sam.hpp"

using gfase::rgb_to_hex;

using BamTools::BamAlignment;
using BamTools::BamReader;


namespace gfase{

class CigarSummary{
public:
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
        path query_path,
        size_t n_threads
);


}


#endif //GFASE_PHASE_ASSIGN_HPP