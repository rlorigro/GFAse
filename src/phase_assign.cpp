#include "phase_assign.hpp"

#include <iomanip>
#include <map>

using ghc::filesystem::create_directories;
using std::setprecision;
using std::min;
using std::map;

namespace gfase{


void CigarSummary::update(char cigar_operation, uint32_t length, uint32_t max_indel) {
    if (cigar_operation == 'M'){
        throw runtime_error("ERROR: M operation must be explicit as '=' or 'X', cannot use ambiguous cigar");
    }
    else if (cigar_operation == '='){
        n_matches += length;

        ref_length += length;
        query_length += length;
    }
    else if (cigar_operation == 'X'){
        n_mismatches += length;

        ref_length += length;
        query_length += length;
    }
    else if (cigar_operation == 'I'){
        n_inserts += min(max_indel, length);

        query_length += min(max_indel, length);
    }
    else if (cigar_operation == 'D'){
        n_deletes += min(max_indel, length);

        ref_length += min(max_indel, length);
    }
}


CigarSummary::CigarSummary():
        n_matches(0),
        n_mismatches(0),
        n_deletes(0),
        n_inserts(0),
        ref_length(0),
        query_length(0)
{}


double CigarSummary::get_identity() {
    return double(n_matches) / (double(n_matches + n_mismatches + n_deletes + n_inserts) + double(1e-12));
}


ostream& operator<<(ostream& o, const CigarSummary& c){
    o << "=:" << '\t' << c.n_matches << '\n';
    o << "X:" << '\t' << c.n_mismatches << '\n';
    o << "I:" << '\t' << c.n_inserts << '\n';
    o << "D:" << '\t' << c.n_deletes << '\n';
    o << "ref:" << '\t' << c.ref_length << '\n';
    o << "query:" << '\t' << c.query_length;

    return o;
}


void parse_bam_cigars(path bam_path, unordered_map<string,CigarSummary>& cigar_summaries, const string& required_prefix){
    BamReader reader;

    if (!reader.Open(bam_path) ) {
        throw std::runtime_error("ERROR: could not read BAM file: " + bam_path.string());
    }

    auto reference_data = reader.GetReferenceData();

    BamAlignment e;

    size_t l = 0;
    while (reader.GetNextAlignment(e) ) {
        if (e.IsPrimaryAlignment() and e.IsMapped()){
            bool valid_prefix = true;
            if (not required_prefix.empty()){
                for (size_t i=0; i<required_prefix.size(); i++){
                    if (e.Name[i] != required_prefix[i]){
                        valid_prefix = false;
                        break;
                    }
                }
            }

            if (not valid_prefix){
                continue;
            }

            if (not e.IsSupplementary()){
                auto& ref_name = reference_data.at(e.RefID).RefName;
                cigar_summaries[e.Name].primary_ref = ref_name;
            }

            for (auto& c: e.CigarData){
                cigar_summaries[e.Name].update(c.Type, c.Length, 20);
            }
        }
        l++;
    }
}


void assign_phases(
        path output_dir,
        path pat_ref_path,
        path mat_ref_path,
        path query_path,
        string required_prefix,
        bool extract_fasta,
        size_t n_threads
){
    create_directories(output_dir);

    // Do some cursory checks on input files
    for (auto& p: {pat_ref_path, mat_ref_path, query_path}){
        auto e = p.extension();
        if (not (e == ".fa" or e == ".fna" or e == ".fasta")){
            throw runtime_error("ERROR: input file does not have extension consistent with FASTA format: " + p.string());
        }

        if (not exists(p)){
            throw runtime_error("ERROR: input file can not be found: " + p.string());
        }
        else{
            ifstream f(p);

            if (not f.good() or not f.is_open()){
                throw runtime_error("ERROR: cannot read file: " + p.string());
            }
        }
    }

    unordered_map<string,CigarSummary> pat_cigar_summaries;
    unordered_map<string,CigarSummary> mat_cigar_summaries;

    map<string,size_t> query_lengths;
    get_query_lengths_from_fasta(query_path, query_lengths);

    // Optionally remove entries that don't contain a prefix specified by user
    vector<string> to_be_deleted;
    for (auto& [name,l]: query_lengths){
        bool valid_prefix = true;
        if (not required_prefix.empty()){
            for (size_t i=0; i<required_prefix.size(); i++){
                if (name[i] != required_prefix[i]){
                    valid_prefix = false;
                    break;
                }
            }
        }

        if (not valid_prefix){
            to_be_deleted.emplace_back(name);
        }
    }

    for (auto& name: to_be_deleted){
        query_lengths.erase(name);
    }

    to_be_deleted.clear();

    // Initialize the cigar summaries with all the input reads, so that there are no missing maps later in the event
    // that it was unmapped in one or both of the reference haplotypes
    for (auto& [name, length]: query_lengths){
        cerr << name << ' ' << length << '\n';
        const CigarSummary c;
        pat_cigar_summaries.emplace(name, c);
        mat_cigar_summaries.emplace(name, c);
    }

    auto query_vs_pat_sam = align(output_dir, pat_ref_path, query_path, n_threads);
    auto query_vs_mat_sam = align(output_dir, mat_ref_path, query_path, n_threads);

    auto query_vs_pat_bam = sam_to_sorted_bam(query_vs_pat_sam, n_threads);
    auto query_vs_mat_bam = sam_to_sorted_bam(query_vs_mat_sam, n_threads);

    parse_bam_cigars(query_vs_pat_bam, pat_cigar_summaries, required_prefix);
    parse_bam_cigars(query_vs_mat_bam, mat_cigar_summaries, required_prefix);

    path output_path = output_dir / "phase_assignments.csv";
    ofstream output_file(output_path);

    auto c0 = rgb_to_hex(0.867,0.133,0.133);
    auto c1 = rgb_to_hex(0.078,0.518,0.518);

    string header = {"name,length,phase,primary_ref,mat_identity,pat_identity,color"};
    output_file << header << '\n';

    unordered_map<string,bool> phased_contigs;

    for (const auto& [name, length]: query_lengths){
        cerr << name << '\n';
        cerr << "pat" << '\n';
        cerr << pat_cigar_summaries.at(name) << '\n';
        cerr << "identity:" << '\t' << std::fixed << std::setprecision(6) << pat_cigar_summaries.at(name).get_identity() << '\n';
        cerr << "mat" << '\n';
        cerr << mat_cigar_summaries.at(name) << '\n';
        cerr << "identity:" << '\t' << std::fixed << std::setprecision(6) << mat_cigar_summaries.at(name).get_identity() << '\n';
        cerr << '\n';

        auto& mat_result = mat_cigar_summaries.at(name);
        auto& pat_result = pat_cigar_summaries.at(name);

        double mat_identity = mat_result.get_identity();
        double pat_identity = pat_result.get_identity();

        bool phase;
        string color;
        string primary_ref;

        if (mat_identity > pat_identity) {
            phase = 1;
            color = "#" + c1;
            primary_ref = mat_result.primary_ref;
        }
        else{
            phase = 0;
            color = "#" + c0;
            primary_ref = pat_result.primary_ref;
        }

        phased_contigs.emplace(name, phase);

        vector<string> s = {
                name,
                to_string(length),
                to_string(phase),
                primary_ref,
                to_string(mat_identity),
                to_string(pat_identity),
                color};

        output_file << join(s, ',') << '\n';

    }



}

}
