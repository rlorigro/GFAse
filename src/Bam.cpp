#include "Bam.hpp"

#include <stdexcept>
#include <iostream>
#include <vector>

using std::runtime_error;
using std::vector;
using std::cerr;
using std::min;


namespace gfase{


Bam::Bam(path bam_path):
    bam_path(bam_path),
    bam_file(nullptr),
    bam_iterator(nullptr)
{
    if ((bam_file = hts_open(bam_path.string().c_str(), "r")) == 0) {
        throw runtime_error("ERROR: Cannot open bam file: " + bam_path.string());
    }

    // bam header
    if ((bam_header = sam_hdr_read(bam_file)) == 0){
        throw runtime_error("ERROR: Cannot open header for bam file: " + bam_path.string() + "\n");
    }

    alignment = bam_init1();
}


void Bam::for_alignment_in_bam(const function<void(const string& ref_name, const string& query_name, uint8_t map_quality, uint16_t flag)>& f){
    while (sam_read1(bam_file, bam_header, alignment) >= 0){
        string query_name = bam_get_qname(alignment);

        string ref_name;

        // Ref name field might be empty if read is unmapped, in which case the target (aka ref) id might not be in range
        if (alignment->core.tid < bam_header->n_targets and alignment->core.tid > 0) {
            ref_name = bam_header->target_name[alignment->core.tid];
        }

        f(ref_name, query_name, alignment->core.qual, alignment->core.flag);
    }
}


void Bam::for_alignment_in_bam(bool get_cigar, const function<void(SamElement& a)>& f){
    while (sam_read1(bam_file, bam_header, alignment) >= 0){
        SamElement e;
        e.query_name = bam_get_qname(alignment);

        // Ref name field might be empty if read is unmapped, in which case the target (aka ref) id might not be in range
        if (alignment->core.tid < bam_header->n_targets and alignment->core.tid > -1) {
            e.ref_name = bam_header->target_name[alignment->core.tid];
        }

        e.mapq = alignment->core.qual;
        e.flag = alignment->core.flag;

        if (get_cigar) {
            auto n_cigar = alignment->core.n_cigar;
            auto cigar_ptr = bam_get_cigar(alignment);
            e.cigars.assign(cigar_ptr, cigar_ptr + n_cigar);
        }

        f(e);
    }
}


void Bam::for_alignment_in_bam(const function<void(FullAlignmentBlock& a)>& f){
    while (sam_read1(bam_file, bam_header, alignment) >= 0){
        FullAlignmentBlock e;

        // Ref name field might be empty if read is unmapped, in which case the target (aka ref) id might not be in range
        if (alignment->core.tid < bam_header->n_targets and alignment->core.tid > -1) {
            e.ref_name = bam_header->target_name[alignment->core.tid];
        }

        e.query_name = bam_get_qname(alignment);
        e.query_length = alignment->core.l_qseq;
        e.mapq = alignment->core.qual;
        e.flag = alignment->core.flag;
        e.ref_start = alignment->core.pos;
        e.is_reverse = bam_is_rev(alignment);

        auto n_cigar = alignment->core.n_cigar;
        auto cigar_ptr = bam_get_cigar(alignment);
        e.cigars.assign(cigar_ptr, cigar_ptr + n_cigar);

        int32_t start_clip = 0;
        int32_t end_clip = 0;
        int32_t query_length;
        e.n_matches = 0;
        e.n_mismatches = 0;
        e.n_inserts = 0;
        e.n_deletes = 0;
        e.n_n = 0;

        size_t i = 0;

        for (auto& c: e.cigars){
            char type = bam_cigar_opchr(c);
            auto length = bam_cigar_oplen(c);

            if (type == 'M' or type == '='){
                e.n_matches += length;
            }
            else if (type == 'X' ){
                e.n_mismatches += length;
            }
            else if (type == 'I'){
                e.n_inserts += length;
            }
            else if (type == 'D'){
                e.n_deletes += length;
            }
            else if (type == 'N'){
                e.n_n += length;
            }
            else if (type == 'S' or type == 'H'){
                if (i == 0){
                    start_clip = int32_t(length);
                }
                else{
                    end_clip = int32_t(length);
                }
            }

            i++;
        }

        query_length = int32_t(e.n_matches) + int32_t(e.n_mismatches) + int32_t(e.n_inserts) + int32_t(start_clip) + int32_t(end_clip);
        e.ref_stop = int32_t(e.ref_start) + int32_t(e.n_matches) + int32_t(e.n_mismatches) + int32_t(e.n_deletes);

        if (e.is_reverse){
            e.query_start = int32_t(end_clip);
            e.query_stop = int32_t(query_length) - int32_t(start_clip);
        }
        else{
            e.query_start = int32_t(start_clip);
            e.query_stop = int32_t(query_length) - int32_t(end_clip);
        }

        f(e);
    }
}


void Bam::for_ref_in_header(const function<void(const string& ref_name, uint32_t length)>& f) const{
    for (size_t i=0; i < bam_header->n_targets; i++){
        const auto& name = bam_header->target_name[i];
        auto length = bam_header->target_len[i];

        f(name, length);
    }
}


Bam::~Bam() {
    hts_close(bam_file);
    bam_hdr_destroy(bam_header);
    bam_destroy1(alignment);
//    hts_idx_destroy(bam_index);
    hts_itr_destroy(bam_iterator);
}


bool Bam::is_first_mate(uint16_t flag){
    return (uint16_t(flag) >> 6) & uint16_t(1);
}


bool Bam::is_second_mate(uint16_t flag){
    return (uint16_t(flag) >> 7) & uint16_t(1);
}


bool Bam::is_not_primary(uint16_t flag){
    return (uint16_t(flag) >> 8) & uint16_t(1);
}


bool Bam::is_primary(uint16_t flag){
    return (not is_not_primary(flag));
}


bool Bam::is_supplementary(uint16_t flag){
    return (uint16_t(flag) >> 11) & uint16_t(1);
}


template <class T> void update_contact_map(
        const vector<T>& alignments,
        MultiContactGraph& contact_graph,
        IncrementalIdMap<string>& id_map){

    // Iterate one triangle of the all-by-all matrix, adding up mapqs for reads on both end of the pair
    for (size_t i=0; i<alignments.size(); i++){
        auto& a = alignments[i];
        auto ref_id_a = int32_t(id_map.try_insert(a.ref_name));
        contact_graph.try_insert_node(ref_id_a, 0);

        contact_graph.increment_coverage(ref_id_a, 1);

        for (size_t j=i+1; j<alignments.size(); j++) {
            auto& b = alignments[j];
            auto ref_id_b = int32_t(id_map.try_insert(b.ref_name));
            contact_graph.try_insert_node(ref_id_b, 0);
            contact_graph.try_insert_edge(ref_id_a, ref_id_b);

            // Increment by the mapq
            if (a.mapq <= 40) {
                contact_graph.increment_edge_weight(ref_id_a, ref_id_b, a.mapq);
            }
            else{
                contact_graph.increment_edge_weight(ref_id_a, ref_id_b, 40);
            }
        }
    }
}


void parse_unpaired_bam_file(
        path bam_path,
        MultiContactGraph& contact_graph,
        IncrementalIdMap<string>& id_map,
        string required_prefix,
        int8_t min_mapq){

    Bam reader(bam_path);

    size_t l = 0;
    string prev_query_name = "";
    vector<SamElement> alignments;

    reader.for_alignment_in_bam(false, [&](const SamElement& a){
        if (l == 0){
            prev_query_name = a.query_name;
        }

        if (prev_query_name != a.query_name){
            update_contact_map(alignments, contact_graph, id_map);
            alignments.clear();
        }

        // No information about reference contig, this alignment is unusable
        if (a.ref_name.empty()){
            return;
        }

        // Optionally filter by the contig names. E.g. "PR" in shasta
        bool valid_prefix = true;
        if (not required_prefix.empty()){
            for (size_t i=0; i<required_prefix.size(); i++){
                if (a.ref_name[i] != required_prefix[i]){
                    valid_prefix = false;
                    break;
                }
            }
        }

        if (valid_prefix) {
            // Only allow reads with mapq > min_mapq and not secondary
            if (a.mapq >= min_mapq and a.is_primary()) {
                alignments.emplace_back(a);
            }
        }

        l++;
        prev_query_name = a.query_name;
    });
}



}
