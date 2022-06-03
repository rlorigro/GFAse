#include "BubbleGraph.hpp"
#include "IncrementalIdMap.hpp"
#include "Filesystem.hpp"
#include "BamReader.h"
#include "CLI11.hpp"
#include "Sam.hpp"

using BamTools::BamReader;
using BamTools::BamAlignment;

using gfase::BubbleGraph;
using gfase::IncrementalIdMap;
using gfase::unpaired_mappings_t;
using gfase::unpaired_mappings_t;
using gfase::SamElement;

using ghc::filesystem::path;
using CLI::App;

#include <iostream>
#include <fstream>
#include <map>

using std::ofstream;
using std::cerr;
using std::min;
using std::map;

using weighted_contact_map_t = sparse_hash_map <int32_t, sparse_hash_map<int32_t, map <int8_t, int32_t> > >;


void parse_unpaired_bam_file(
        path bam_path,
        unpaired_mappings_t& mappings,
        IncrementalIdMap<string>& id_map,
        string required_prefix,
        int8_t min_mapq){

    BamReader reader;

    if (!reader.Open(bam_path) ) {
        throw std::runtime_error("ERROR: could not read BAM file: " + bam_path.string());
    }

    auto reference_data = reader.GetReferenceData();

    BamAlignment a;
    size_t l = 0;

    while (reader.GetNextAlignment(a) ) {
//        cerr << "e.Name" << ' ' << a.Name << '\n';
//        cerr << "e.RefID" << ' ' << a.RefID << '\n';
//        cerr << "int32_t(l)" << ' ' << int32_t(l) << '\n';
//        cerr << "int(e.AlignmentFlag)" << ' ' << int(a.AlignmentFlag) << '\n';
//        cerr << "int(e.IsPrimaryAlignment)" << ' ' << int(a.IsPrimaryAlignment()) << '\n';
//        cerr << "int(e.IsSupplementary)" << ' ' << int(a.IsSupplementary()) << '\n';
//        cerr << "int(e.MapQuality)" << ' ' << int(a.MapQuality) << '\n';
//        cerr << '\n';

        // No information about reference contig, this alignment is unusable
        if (a.RefID < 0 or a.RefID > reference_data.size()){
            continue;
        }

        auto& ref_name = reference_data.at(a.RefID).RefName;

        bool valid_prefix = true;
        if (not required_prefix.empty()){
            for (size_t i=0; i<required_prefix.size(); i++){
                if (ref_name[i] != required_prefix[i]){
                    valid_prefix = false;
                    break;
                }
            }
        }

        if (valid_prefix) {
            id_map.try_insert(ref_name);

            if (a.MapQuality >= min_mapq and a.IsPrimaryAlignment()) {
                SamElement s(a.Name, ref_name, int32_t(l), int16_t(a.AlignmentFlag), int8_t(a.MapQuality));
                mappings[a.Name].emplace(s);
            }
        }

        l++;
    }
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


void generate_contact_map_from_mappings(
        const unpaired_mappings_t& mappings,
        const IncrementalIdMap<string>& id_map,
        weighted_contact_map_t& contact_map
){
    // All-vs-all, assuming reference names are unique in first vs second mates
    for (const auto& [name,elements]: mappings) {
        size_t i = 0;

        for (const auto& e: elements){
            size_t j = 0;

            for (const auto& e2: elements){
                if (j > i){
                    if (e.ref_name != e2.ref_name) {
                        auto id = int32_t(id_map.get_id(e.ref_name));
                        auto id2 = int32_t(id_map.get_id(e2.ref_name));

                        auto min_mapq = min(e.mapq, e2.mapq);

                        contact_map[id][id2][min_mapq]++;
                        contact_map[id2][id][min_mapq]++;
                    }
                }
                j++;
            }
            i++;
        }
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
            output_file << id_map.get_name(id) << ',' << id_map.get_name(id2) << ',' << '{';

            auto max_key = map3.rbegin()->first;

            for (const auto& [q,count]: map3) {
                 output_file << int(q) << ':' << count;

                 if (q != max_key){
                    output_file << ',';
                 }
            }
            output_file << '}' << '\n';
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

    // TODO: Load GFA first, and find a method for identifying haplotypic bubbles
    // TODO: in other words, stop relying on shasta conventions

    // Mappings grouped by their read name (to simplify paired-end parsing)
    unpaired_mappings_t mappings;

    // Datastructures to represent linkages from hiC
    weighted_contact_map_t contact_map;
    vector <vector <int32_t> > adjacency;

    cerr << "Loading alignments..." << '\n';

    if (sam_path.extension() == ".bam"){
        parse_unpaired_bam_file(sam_path, mappings, id_map, required_prefix, min_mapq);
    }
    else{
        throw runtime_error("ERROR: unrecognized extension for SAM/BAM input file: " + sam_path.extension().string());
    }

    remove_singleton_reads(mappings);

    cerr << "Generating contact map..." << '\n';

    // Build the contact map by iterating sets of alignments and creating edges in an all-by-all fashion within sets
    generate_contact_map_from_mappings(mappings, id_map, contact_map);

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

