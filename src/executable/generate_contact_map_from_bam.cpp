#include "BubbleGraph.hpp"
#include "IncrementalIdMap.hpp"
#include "Filesystem.hpp"
#include "CLI11.hpp"
#include "Sam.hpp"
#include "Bam.hpp"

using gfase::BubbleGraph;
using gfase::IncrementalIdMap;
using gfase::unpaired_mappings_t;
using gfase::unpaired_mappings_t;
using gfase::SamElement;
using gfase::Bam;

using ghc::filesystem::path;
using CLI::App;

#include <iostream>
#include <fstream>
#include <map>

using std::ofstream;
using std::cerr;
using std::min;
using std::map;

using weighted_contact_map_t = sparse_hash_map <int32_t, sparse_hash_map<int32_t, map <uint8_t, int32_t> > >;
using mappings_per_read_t = sparse_hash_map <string, map <size_t, map <uint8_t, int64_t> > >;


void update_contact_map(
        vector<SamElement>& alignments,
        weighted_contact_map_t& contact_map,
        IncrementalIdMap<string>& id_map){

    // Iterate one triangle of the all-by-all matrix, adding up mapqs for reads on both end of the pair
    for (size_t i=0; i<alignments.size(); i++){
//        cerr << alignments[i] << '\n';

        for (size_t j=i+1; j<alignments.size(); j++) {
            auto& a = alignments[i];
            auto& b = alignments[j];

            auto ref_id_a = id_map.try_insert(a.ref_name);
            auto ref_id_b = id_map.try_insert(b.ref_name);

            // TODO: split left and right mapq?
            contact_map[ref_id_a][ref_id_b][min(a.mapq,b.mapq)]++;
            contact_map[ref_id_b][ref_id_a][min(a.mapq,b.mapq)]++;
        }
    }
}


void parse_unpaired_bam_file(
        path bam_path,
        weighted_contact_map_t& contact_map,
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
            if (alignments.size() > 1){
                update_contact_map(alignments, contact_map, id_map);
            }
            alignments.clear();
        }

        // No information about reference contig, this alignment is unusable
        if (a.ref_name.empty()){
            return;
        }

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
            if (a.mapq >= min_mapq and a.is_primary()) {
                alignments.emplace_back(a);
            }
        }

        l++;
        prev_query_name = a.query_name;
    });
}


void remove_singleton_reads(unpaired_mappings_t& mappings){
    vector<string> to_be_deleted;
    for (auto& [name,elements]: mappings){

        // Remove all the entries where only one mapping exists for that read
        if (elements.size() < 2){
            to_be_deleted.emplace_back(name);
        }
    }

    for (auto& item: to_be_deleted){
        mappings.erase(item);
    }
}


void write_contact_map(
        path output_path,
        const weighted_contact_map_t& contact_map,
        const IncrementalIdMap<string>& id_map){
    ofstream output_file(output_path);

    if (not output_file.is_open() or not output_file.good()){
        throw std::runtime_error("ERROR: could not write to file: " + output_path.string());
    }

    for (const auto& [id,map2]: contact_map){
        for (const auto& [id2,map3]: map2){
            output_file << id_map.get_name(id) << ',' << id_map.get_name(id2) << ',';

            auto max_key = map3.rbegin()->first;

            for (const auto& [q,count]: map3) {
                 output_file << int(q) << ':' << count;

                 if (q != max_key){
                    output_file << ' ';
                 }
            }
            output_file << '\n';
        }
    }
}


void generate_contact_map_from_bam(path output_dir, path sam_path, path gfa_path, string required_prefix, int8_t min_mapq, size_t n_threads){
    if (exists(output_dir)){
        throw runtime_error("ERROR: output directory exists already");
    }
    else {
        create_directories(output_dir);
    }

    // Id-to-name bimap for reference contigs
    IncrementalIdMap<string> id_map(false);

    // Datastructures to represent linkages from hiC
    weighted_contact_map_t contact_map;
    vector <vector <int32_t> > adjacency;

    cerr << "Loading alignments as contact map..." << '\n';

    if (sam_path.extension() == ".bam"){
        parse_unpaired_bam_file(sam_path, contact_map, id_map, required_prefix, min_mapq);
    }
    else{
        throw runtime_error("ERROR: unrecognized extension for SAM/BAM input file: " + sam_path.extension().string());
    }

    path output_path = output_dir / "contacts.csv";
    write_contact_map(output_path, contact_map, id_map);
}


int main (int argc, char* argv[]){
    path sam_path;
    path gfa_path;
    path output_dir;
    string required_prefix;
    int8_t min_mapq = 0;
    size_t n_threads = 1;

    CLI::App app{"App description"};

    app.add_option(
            "-i,--input",
            sam_path,
            "Path to SAM or BAM containing filtered, paired HiC reads")
            ->required();

    app.add_option(
            "-g,--gfa",
            gfa_path,
            "Path to GFA containing assembly graph to be phased");

    app.add_option(
            "-o,--output_dir",
            output_dir,
            "Path to (nonexistent) directory where output will be stored")
            ->required();

    app.add_option(
            "-p,--prefix",
            required_prefix,
            "Prefix required in ref name for mapping to be counted");

    app.add_option(
            "-m,--min_mapq",
            min_mapq,
            "Minimum required mapq value for mapping to be counted");

    app.add_option(
            "-t,--threads",
            n_threads,
            "Maximum number of threads to use");

    CLI11_PARSE(app, argc, argv);

    generate_contact_map_from_bam(output_dir, sam_path, gfa_path, required_prefix, min_mapq, n_threads);

    return 0;
}

