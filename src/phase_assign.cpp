#include "phase_assign.hpp"
#include "misc.hpp"

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
            cerr << "BAM NAME " << e.Name << '\n';

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


void bin_fasta_sequence(string& name, string& sequence, int bin, ofstream& file0, ofstream& file1){
    cerr << name << ' ' << sequence.substr(0,min(sequence.size(),size_t(10))) << '\n';

    if (not sequence.empty()){
        if (bin == 0){
            file0 << '>' << name << '\n';
            file0 << sequence << '\n';
        }
        else if (bin == 1){
            file1 << '>' << name << '\n';
            file1 << sequence << '\n';
        }
    }
}


void bin_fasta_sequences(path input_fasta_path, path output_pat_fasta_path, path output_mat_fasta_path, const array <set <string>, 2>& phased_contigs){
    ifstream input_file(input_fasta_path);
    ofstream file0(output_pat_fasta_path);
    ofstream file1(output_mat_fasta_path);

    string name;
    string sequence;

    // 0=pat, 1=mat, 2=both
    int bin = 2;
    char c;

    while (input_file.get(c)){
        if (c == '>'){
            bin_fasta_sequence(name, sequence, bin, file0, file1);

            sequence.clear();
            name.clear();

            while (input_file.get(c)){
                if (std::isspace(c)){
                    // Check if this segment is phased, set flag accordingly
                    auto count_0 = phased_contigs[0].count(name);
                    auto count_1 = phased_contigs[1].count(name);

                    if ((count_0 > 0) and (count_1 == 0)){
                        bin = 0;
                    }
                    else if ((count_0 == 0) and (count_1 > 0)){
                        bin = 1;
                    }
                    else if ((count_0 > 0) and (count_1 > 0)){
                        throw runtime_error("ERROR: contig assigned phase 0 and phase 1: " + name);
                    }
                    else{
                        bin = 2;
                    }
                    break;
                }
                else {
                    name += c;
                }
            }

            // Skip rest of the header line which may contain space-separated data
            input_file.ignore(std::numeric_limits<streamsize>::max(), '\n');
        }
        else if (c != '\n'){
            sequence += c;
        }
    }

    bin_fasta_sequence(name, sequence, bin, file0, file1);
}


void assign_phases(
        path output_dir,
        path pat_ref_path,
        path mat_ref_path,
        path query_path,
        string required_prefix,
        size_t n_threads,
        array <set <string>, 2>& phased_contigs,
        map<string,size_t>& query_lengths,
        bool extract_fasta
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
    for (const auto& [name, length]: query_lengths){
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

        int phase = 2;
        string color;
        string primary_ref;

        if (mat_identity != 0 and pat_identity != 0){
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

            phased_contigs[phase].emplace(name);
        }
        else{
            color = "#708090";
        }

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

    if (extract_fasta) {
        path pat_output_fasta_path = output_dir / "pat_sequences.fasta";
        path mat_output_fasta_path = output_dir / "mat_sequences.fasta";

        bin_fasta_sequences(query_path, pat_output_fasta_path, mat_output_fasta_path, phased_contigs);
    }
}


/*
    #####  Color Palette by Paletton.com
    #####  Palette URL: http://paletton.com/#uid=7000u0kr6rGgDzJlztwt9l+yoh7

    *** Primary color:

       shade 0 = #DD2222 = rgb(221, 34, 34) = rgba(221, 34, 34,1) = rgb0(0.867,0.133,0.133)

    *** Secondary color (1):

       shade 0 = #DD7622 = rgb(221,118, 34) = rgba(221,118, 34,1) = rgb0(0.867,0.463,0.133)

    *** Secondary color (2):

       shade 0 = #148484 = rgb( 20,132,132) = rgba( 20,132,132,1) = rgb0(0.078,0.518,0.518)

    *** Complement color:

       shade 0 = #1BB01B = rgb( 27,176, 27) = rgba( 27,176, 27,1) = rgb0(0.106,0.69,0.106)

    #####  Generated by Paletton.com (c) 2002-2014
 */
void write_bandage_csv(
        const vector<string>& intersection_00,
        const vector<string>& intersection_11,
        const vector<string>& intersection_01,
        const vector<string>& intersection_10,
        path output_path){

    ofstream file(output_path);

    file << "Name" << ',' << "Set" << ',' << "Color" << '\n';
    for (const auto& item: intersection_00){
        file << item << ',' << "00" << ',' << '#' << rgb_to_hex(0.867,0.133,0.133) << '\n';     // Red
    }
    for (const auto& item: intersection_11){
        file << item << ',' << "11" << ',' << '#' << rgb_to_hex(0.867,0.463,0.133) << '\n';     // Orange
    }
    for (const auto& item: intersection_01){
        file << item << ',' << "01" << ',' << '#' << rgb_to_hex(0.078,0.518,0.518) << '\n';     // Blue
    }
    for (const auto& item: intersection_10){
        file << item << ',' << "10" << ',' << '#' << rgb_to_hex(0.106,0.69,0.106) << '\n';      // Green
    }
}


size_t get_total_length_of_phase_set(const map<string,size_t>& contig_lengths, const vector<string>& s){
    size_t l = 0;

    for (const auto& item: s){
        l += contig_lengths.at(item);
    }

    return l;
}


void evaluate_phasing(
        path output_dir,
        path contact_phase_csv,
        path pat_ref_path,
        path mat_ref_path,
        path query_path,
        string required_prefix,
        size_t n_threads,
        bool extract_fasta
){
    map<string,size_t> query_lengths;

    array <set <string>, 2> inferred_phases;
    array <set <string>, 2> true_phases;

    for_entry_in_csv(contact_phase_csv, [&](const vector<string>& tokens, size_t line){
        // Skip header
        if (line == 0){
            return;
        }

        bool phase;

        try {
            phase = stoi(tokens[1]);
        }
        catch (exception& e){
            cerr << e.what() << '\n';
            throw std::runtime_error("ERROR: tokens[1] cannot be parsed: " + tokens[0] + ' ' + tokens[1] + '\n');
        }

        auto name = tokens[0];


        inferred_phases[phase].emplace(name);
    });

    assign_phases(
            output_dir,
            pat_ref_path,
            mat_ref_path,
            query_path,
            required_prefix,
            n_threads,
            true_phases,
            query_lengths,
            extract_fasta);

    vector<string> intersection_00;
    vector<string> intersection_11;
    vector<string> intersection_01;
    vector<string> intersection_10;

    set_intersection(true_phases[0].begin(), true_phases[0].end(),
                     inferred_phases[0].begin(), inferred_phases[0].end(),
                     std::back_inserter(intersection_00));

    set_intersection(true_phases[1].begin(), true_phases[1].end(),
                     inferred_phases[1].begin(), inferred_phases[1].end(),
                     std::back_inserter(intersection_11));

    set_intersection(true_phases[0].begin(), true_phases[0].end(),
                     inferred_phases[1].begin(), inferred_phases[1].end(),
                     std::back_inserter(intersection_01));

    set_intersection(true_phases[1].begin(), true_phases[1].end(),
                     inferred_phases[0].begin(), inferred_phases[0].end(),
                     std::back_inserter(intersection_10));

    cerr << "Intersection" << '\n';
    cerr << '\t' << "00" << ' ' << intersection_00.size() << ' ' << get_total_length_of_phase_set(query_lengths, intersection_00) << '\n';
    cerr << '\t' << "11" << ' ' << intersection_11.size() << ' ' << get_total_length_of_phase_set(query_lengths, intersection_11) << '\n';
    cerr << '\t' << "01" << ' ' << intersection_01.size() << ' ' << get_total_length_of_phase_set(query_lengths, intersection_01) << '\n';
    cerr << '\t' << "10" << ' ' << intersection_10.size() << ' ' << get_total_length_of_phase_set(query_lengths, intersection_10) << '\n';

    write_bandage_csv(intersection_00, intersection_11, intersection_01, intersection_10, output_dir / "phase_intersections.csv");

}


}
