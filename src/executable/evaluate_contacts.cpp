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
using gfase::edge;
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
        Bam& reader,
        const function<void(const FullAlignmentChain& read_alignments)>& f
        ){
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


class Histogram{
public:
    vector<int32_t> frequencies;
    size_t max_value;
    size_t bin_size;

    Histogram(size_t max_value, size_t bin_size);
    void for_each_item(const function<void(size_t value, int32_t frequency)>& f) const;
    void get_reverse_cdf(vector<int32_t>& cdf) const;
    void update(size_t value);
    void update(size_t value, int32_t count);
    void write_to_csv(path output_path) const;
};


Histogram::Histogram(size_t max_value, size_t bin_size):
    frequencies(max_value/bin_size + 1),
    max_value(max_value),
    bin_size(bin_size)
{}


void Histogram::update(size_t value){
    if (value <= max_value){
        frequencies[value/bin_size]++;
    }
}


void Histogram::update(size_t value, int32_t count){
    if (value <= max_value){
        frequencies[value/bin_size]+= count;
    }
}


void Histogram::for_each_item(const function<void(size_t value, int32_t frequency)>& f) const{
    for (size_t i=0; i<frequencies.size(); i++){
        f(i*bin_size, frequencies[i]);
    }
}


void Histogram::write_to_csv(path output_path) const{
    ofstream file(output_path);

    if (not (file.is_open() and file.good())){
        throw runtime_error("ERROR: could not write to file: " + output_path.string());
    }

    for (size_t i=0; i<frequencies.size(); i++){
        if (frequencies[i] == 0){
            continue;
        }
        file << i*bin_size << ',' << frequencies[i] << '\n';
    }
}



void Histogram::get_reverse_cdf(vector<int32_t>& cdf) const{
    cdf.clear();
    cdf.resize(frequencies.size());

    int64_t i = int64_t(frequencies.size()) - 1;
    int32_t prev = 0;
    for (auto it = frequencies.rbegin(); it<frequencies.rend(); it++){
        cdf[i] = *it + prev;
        prev = cdf[i];
        i--;
    }
}


class StratifiedContactMap{
public:
    unordered_map <pair <int32_t, int32_t>, Histogram> contacts;
    unordered_map <int32_t, int8_t> partitions;

    StratifiedContactMap()=default;
    void update(int32_t a, int32_t b, uint8_t mapq);
    void set_partition(int32_t id, int8_t partition);
    void for_each_edge(const function<void(const pair<int32_t,int32_t>& edge, const Histogram& weights)>& f) const;
    int8_t get_partition(int32_t id) const;
};


void StratifiedContactMap::update(int32_t a, int32_t b, uint8_t mapq){
    auto e = edge(a,b);

    auto result = contacts.find(e);

    // Make a new histogram if this edge doesn't exist yet
    if (result == contacts.end()){
        result = contacts.emplace(e,Histogram(60,1)).first;
    }

    // Update the (now guaranteed to exist) histogram
    result->second.update(mapq);
}


void StratifiedContactMap::set_partition(int32_t id, int8_t partition){
    partitions[id] = partition;
}


int8_t StratifiedContactMap::get_partition(int32_t id) const{
    return partitions.at(id);
}


void StratifiedContactMap::for_each_edge(const function<void(const pair<int32_t,int32_t>& edge, const Histogram& weights)>& f) const{
    for (auto& [e,w]: contacts){
        f(e,w);
    }
}


void evaluate_contacts(path bam_path, path phase_csv, path output_dir){
    if (exists(output_dir)){
        throw runtime_error("ERROR: output directory exists already");
    }
    else {
        create_directories(output_dir);
    }

    Bam reader(bam_path);
    IncrementalIdMap<string> id_map;
    StratifiedContactMap contact_map;

    reader.for_ref_in_header([&](const string& ref_name, uint32_t length){
        id_map.try_insert(ref_name);
    });

    Histogram subread_lengths(10000000, 50);
    Histogram subread_counts(10000, 1);
    Histogram gap_lengths(1000000000, 500);
    Histogram mapqs(100, 1);

    Histogram n_consistent_contacts(60,1);
    Histogram n_inconsistent_contacts(60,1);
    size_t total_alignments = 0;

    for_each_read_in_bam(reader, [&](const FullAlignmentChain& alignments){
        subread_counts.update(alignments.chain.size());
        total_alignments += alignments.chain.size();

        // Skip processing singleton chains, but note them in the chain length distribution
        if (alignments.chain.size() == 1){
            return;
        }

        // Iterate one triangle of the all-by-all matrix, accumulating stats for alignments and linkages
        for (size_t i=0; i<alignments.chain.size(); i++){
            auto& a = alignments.chain[i];
            auto ref_id_a = int32_t(id_map.get_id(a.ref_name));

            for (size_t j=i+1; j<alignments.chain.size(); j++) {
                auto& b = alignments.chain[j];
                auto ref_id_b = int32_t(id_map.get_id(b.ref_name));

                // Update contact map
                auto min_mapq = min(a.mapq, b.mapq);
                contact_map.update(ref_id_a,ref_id_b,min_mapq);

                // Update gap lengths
                if (ref_id_a == ref_id_b){
                    auto a_middle = (a.ref_stop + a.ref_start)/2;
                    auto b_middle = (b.ref_stop + b.ref_start)/2;
                    auto distance = max(a_middle,b_middle) - min(a_middle,b_middle);

                    gap_lengths.update(distance);
                }
            }

            auto length = a.query_stop - a.query_start;
            subread_lengths.update(length);
            mapqs.update(a.mapq);
        }
    });

    // Read phases from CSV
    for_each_item_in_phase_csv(phase_csv, [&](const string& name, int8_t phase){
        auto id = int32_t(id_map.get_id(name));
        contact_map.set_partition(id, phase);
    });

    // Find which edges are cross-phase or not
    contact_map.for_each_edge([&](const pair<int32_t,int32_t> e, const Histogram& weights){
        auto p_a = contact_map.get_partition(e.first);
        auto p_b = contact_map.get_partition(e.second);

        vector<int32_t> reverse_cdf;
        weights.get_reverse_cdf(reverse_cdf);

        if ((p_a != 0) and (p_b != 0)){
            if (p_a == p_b){
                for (size_t i=0; i<reverse_cdf.size(); i++) {
                    n_consistent_contacts.update(i, reverse_cdf[i]);
                }
            }
            else{
                for (size_t i=0; i<reverse_cdf.size(); i++) {
                    n_inconsistent_contacts.update(i, reverse_cdf[i]);
                }
                auto a_name = id_map.get_name(e.first);
                auto b_name = id_map.get_name(e.second);
            }
        }
    });

    path output_path;

    output_path = output_dir / "subread_lengths.csv";
    subread_lengths.write_to_csv(output_path);
    output_path = output_dir / "subread_counts.csv";
    subread_counts.write_to_csv(output_path);
    output_path = output_dir / "gap_lengths.csv";
    gap_lengths.write_to_csv(output_path);
    output_path = output_dir / "mapq.csv";
    mapqs.write_to_csv(output_path);

    output_path = output_dir / "contacts_summary.csv";
    ofstream contacts_summary_file(output_path);

    if (not (contacts_summary_file.is_open() and contacts_summary_file.good())){
        throw runtime_error("ERROR: could not write to file: " + output_path.string());
    }

    contacts_summary_file << "total_alignments" << ',' << total_alignments << '\n';
    contacts_summary_file << "singletons" << ',' << subread_counts.frequencies[1] << '\n';

    output_path = output_dir / "phasing_summary.csv";
    ofstream phasing_summary_file(output_path);

    if (not (phasing_summary_file.is_open() and phasing_summary_file.good())){
        throw runtime_error("ERROR: could not write to file: " + output_path.string());
    }

    phasing_summary_file << "" << ',';
    n_consistent_contacts.for_each_item([&](size_t value, int32_t frequency){
        phasing_summary_file << value << ',';
    });
    phasing_summary_file << '\n';

    phasing_summary_file << "n_consistent_contacts" << ',';
    n_consistent_contacts.for_each_item([&](size_t value, int32_t frequency){
        phasing_summary_file << frequency << ',';
    });
    phasing_summary_file << '\n';

    phasing_summary_file << "n_inconsistent_contacts" << ',';
    n_inconsistent_contacts.for_each_item([&](size_t value, int32_t frequency){
        phasing_summary_file << frequency << ',';
    });
    phasing_summary_file << '\n';

    phasing_summary_file << "signal_ratio" << ',';
    for (size_t i=0; i<n_consistent_contacts.frequencies.size(); i++){
        auto a = double(n_consistent_contacts.frequencies[i]);
        auto b = double(n_inconsistent_contacts.frequencies[i]);
        phasing_summary_file << a/b << ',';
    }
    phasing_summary_file << '\n';
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

