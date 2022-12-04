#include "MultiContactGraph.hpp"
#include "IncrementalIdMap.hpp"
#include "Filesystem.hpp"
#include "CLI11.hpp"
#include "Bam.hpp"
#include "Sam.hpp"

using gfase::MultiContactGraph;
using gfase::FullAlignmentChain;
using gfase::FullAlignmentBlock;
using gfase::IncrementalIdMap;
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
using std::tuple;
using std::cerr;
using std::map;
using std::min;
using std::max;


// Cross chromosome contacts?
// Cross phase contacts
// Mapq distribution
// Number of contacts
// Alignment length distribution
// Gap size?


void for_each_read_in_bam(
        path bam_path,
        const function<void(const FullAlignmentChain& read_alignments)>& f
        ){

    Bam reader(bam_path);

    size_t l = 0;
    string prev_query_name = "";
    FullAlignmentChain alignments;

    reader.for_alignment_in_bam([&](const FullAlignmentBlock& a){
        if (l == 0){
            prev_query_name = a.query_name;
        }

        if (prev_query_name != a.query_name){
            f(alignments);
            alignments = {};
        }

        // No information about reference contig, this alignment is unusable
        if (a.ref_name.empty()){
            return;
        }

        alignments.chain.emplace_back(a);

        l++;
        prev_query_name = a.query_name;
    });

    // Collect final read's contacts
    if (not alignments.empty()){
        f(alignments);
    }
}


void write_vector_distribution_to_file(const vector<int32_t>& distribution, int32_t bin_size, path output_path){
    ofstream file(output_path);

    if (not (file.is_open() and file.good())){
        throw runtime_error("ERROR: could not write to file: " + output_path.string());
    }

    for (size_t i=0; i<distribution.size(); i++){
        if (distribution[i] == 0){
            continue;
        }
        file << i*bin_size << ',' << distribution[i] << '\n';
    }
}


void for_each_item_in_phase_csv(path csv_path, const function<void(const string& name, int8_t phase)>& f){
    string line;
    string name;
    string phase;

    ifstream file(csv_path);

    if (not (file.is_open() and file.good())){
        throw runtime_error("ERROR: couldn't read file: " + csv_path.string());
    }

    while (getline(file,line)){
        auto pos = line.find_first_of(',');
        name = line.substr(0,pos);
        line = line.substr(pos+1);
        phase = line.substr(0,min(line.size(), line.find_first_of(',')));

        f(name, stoi(phase));
    }
}


void evaluate_contacts(path bam_path, path phase_csv, path output_dir){
    if (exists(output_dir)){
        throw runtime_error("ERROR: output directory exists already");
    }
    else {
        create_directories(output_dir);
    }

    IncrementalIdMap<string> id_map;
    MultiContactGraph contact_graph;
    vector<int32_t> subread_lengths(10000000, 0);
    vector<int32_t> subread_counts(10000, 0);
    vector<int32_t> gap_lengths(10000000, 0);
    vector<int32_t> mapqs(100, 0);

    vector <tuple <string,string,int32_t> > inconsistent_contacts;
    int64_t n_consistent_contacts = 0;
    int64_t n_inconsistent_contacts = 0;
    int64_t total_alignments = 0;

    int32_t bin_size = 100;

    for_each_read_in_bam(bam_path, [&](const FullAlignmentChain& alignments){
        if (alignments.chain.size() < subread_counts.size() - 1){
            subread_counts[alignments.chain.size()]++;
        }

        total_alignments += int64_t(alignments.chain.size());

        // Skip processing singleton chains, but note them in the chain length distribution
        if (alignments.chain.size() == 1){
            return;
        }

        // Iterate one triangle of the all-by-all matrix, accumulating stats for alignments and linkages
        for (size_t i=0; i<alignments.chain.size(); i++){
            auto& a = alignments.chain[i];
            auto ref_id_a = int32_t(id_map.try_insert(a.ref_name));
            contact_graph.try_insert_node(ref_id_a, 0);

            contact_graph.increment_coverage(ref_id_a, 1);

            for (size_t j=i+1; j<alignments.chain.size(); j++) {
                auto& b = alignments.chain[j];
                auto ref_id_b = int32_t(id_map.try_insert(b.ref_name));
                contact_graph.try_insert_node(ref_id_b, 0);
                contact_graph.try_insert_edge(ref_id_a, ref_id_b);
                contact_graph.increment_edge_weight(ref_id_a, ref_id_b, 1);

                if (ref_id_a == ref_id_b){
                    auto a_middle = (a.ref_stop + a.ref_start)/2;
                    auto b_middle = (b.ref_stop + b.ref_start)/2;
                    auto distance = max(a_middle,b_middle) - min(a_middle,b_middle);
                    auto distance_bin = distance/bin_size;

                    if (distance_bin < gap_lengths.size() - 1){
                        gap_lengths[distance_bin]++;
                    }
                }
            }

            auto length = a.query_stop - a.query_start;
            auto length_bin = length/bin_size;

            if (length_bin < subread_lengths.size() - 1){
                subread_lengths[length_bin]++;
            }

            if (a.mapq < mapqs.size() - 1){
                mapqs[a.mapq]++;
            }
        }
    });

    for_each_item_in_phase_csv(phase_csv, [&](const string& name, int8_t phase){
        auto id = int32_t(id_map.try_insert(name));

        if (contact_graph.has_node(id)){
            contact_graph.set_partition(id,phase);
        }
    });

    contact_graph.for_each_edge([&](const pair<int32_t,int32_t> e, int32_t weight){
        auto p_a = contact_graph.get_partition(e.first);
        auto p_b = contact_graph.get_partition(e.second);

        if ((p_a != 0) and (p_b != 0)){
            if (p_a == p_b){
                n_consistent_contacts += weight;
            }
            else{
                n_inconsistent_contacts += weight;
                auto a_name = id_map.get_name(e.first);
                auto b_name = id_map.get_name(e.second);
                inconsistent_contacts.emplace_back(a_name, b_name, weight);
            }
        }
    });

    path output_path;

    cerr << "Writing to file: observed lengths " << '\n';
    output_path = output_dir / "subread_lengths.csv";
    write_vector_distribution_to_file(subread_lengths, bin_size, output_path);

    cerr << "Writing to file: observed gaps " << '\n';
    output_path = output_dir / "subread_counts.csv";
    write_vector_distribution_to_file(subread_counts, 1, output_path);

    cerr << "Writing to file: observed subread counts" << '\n';
    output_path = output_dir / "gap_lengths.csv";
    write_vector_distribution_to_file(gap_lengths, bin_size, output_path);

    cerr << "Writing to file: observed mapq counts" << '\n';
    output_path = output_dir / "mapq.csv";
    write_vector_distribution_to_file(mapqs, 1, output_path);

    cerr << "Writing to file: inconsistent contacts" << '\n';
    output_path = output_dir / "inconsistent_contacts.csv";
    ofstream inconsistent_contacts_file(output_path);

    if (not (inconsistent_contacts_file.is_open() and inconsistent_contacts_file.good())){
        throw runtime_error("ERROR: could not write to file: " + output_path.string());
    }

    for (auto& [a,b,n]: inconsistent_contacts){
        inconsistent_contacts_file << a << ',' << b << ',' << n << '\n';
    }

    cerr << "Writing to file: contacts summary" << '\n';
    output_path = output_dir / "contacts_summary.csv";
    ofstream contacts_summary_file(output_path);

    if (not (contacts_summary_file.is_open() and contacts_summary_file.good())){
        throw runtime_error("ERROR: could not write to file: " + output_path.string());
    }

    contacts_summary_file << "total_alignments" << ',' << total_alignments << '\n';
    contacts_summary_file << "singletons" << ',' << subread_counts[1] << '\n';
    contacts_summary_file << "n_consistent_contacts" << ',' << n_consistent_contacts << '\n';
    contacts_summary_file << "n_inconsistent_contacts" << ',' << n_inconsistent_contacts << '\n';
    contacts_summary_file << "signal_ratio" << ',' << double(n_consistent_contacts)/double(n_inconsistent_contacts) << '\n';
}


int main (int argc, char* argv[]){
    path bam_path;
    path phase_csv;
    path output_dir;

    CLI::App app{"App description"};

    app.add_option(
            "-i,--input",
            bam_path,
            "Path to SAM containing filtered, paired HiC reads")
            ->required();

    app.add_option(
            "-p,--phase_csv",
            phase_csv,
            "Path to CSV indicating phase of contigs")
            ->required();

    app.add_option(
            "-o,--output_dir",
            output_dir,
            "Path to CSV indicating phase of contigs")
            ->required();

    CLI11_PARSE(app, argc, argv);

    evaluate_contacts(bam_path, phase_csv, output_dir);

    return 0;
}

