#include "bamtools/api/BamReader.h"
#include "IncrementalIdMap.hpp"
#include "gfa_to_handle.hpp"
#include "graph_utility.hpp"
#include "BubbleGraph.hpp"
#include "Bipartition.hpp"
#include "hash_graph.hpp"
#include "Filesystem.hpp"
//#include "sparsepp/spp.h"
#include "CLI11.hpp"
#include "Sam.hpp"

using gfase::for_element_in_sam_file;
using gfase::unpaired_mappings_t;
using gfase::paired_mappings_t;
using gfase::contact_map_t;

using gfase::gfa_to_handle_graph;
using gfase::IncrementalIdMap;
using gfase::BubbleGraph;
using gfase::SamElement;

using bdsg::HashGraph;

using BamTools::BamAlignment;
using BamTools::BamReader;

using ghc::filesystem::path;
//using spp::sparse_hash_map;
using CLI::App;

#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <stdexcept>
#include <fstream>
#include <utility>
#include <atomic>
#include <thread>
#include <limits>
#include <bitset>
#include <vector>
#include <mutex>
#include <array>
#include <set>

using std::numeric_limits;
using std::unordered_set;
using std::unordered_map;
using std::runtime_error;
using std::streamsize;
using std::exception;
using std::to_string;
using std::function;
using std::ifstream;
using std::ofstream;
using std::shuffle;
using std::string;
using std::vector;
using std::bitset;
using std::thread;
using std::atomic;
using std::array;
using std::mutex;
using std::pair;
using std::stoi;
using std::cerr;
using std::cref;
using std::ref;
using std::set;


void print_mappings(const unpaired_mappings_t& mappings){
    for (const auto& [name,elements]: mappings){
        cerr << '\n';
        cerr << name << '\n';

        cerr << "Alignments:" << '\n';
        for (auto& e: elements){
            cerr << '\t' << e.ref_name << '\t' << bitset<sizeof(e.flag)*8>(e.flag) << '\t' << "MQ" << int(e.mapq) << ' ' << 'P' << !e.is_not_primary() << ' ' << 'S' << e.is_supplementary() << '\n';
        }
    }
}


void parse_unpaired_sam_file(
        path sam_path,
        unpaired_mappings_t& mappings,
        IncrementalIdMap<string>& id_map,
        string required_prefix,
        int8_t min_mapq){

    for_element_in_sam_file(sam_path, [&](SamElement& e){
        bool valid_prefix = true;
        if (not required_prefix.empty()){
            for (size_t i=0; i<required_prefix.size(); i++){
                if (e.ref_name[i] != required_prefix[i]){
                    valid_prefix = false;
                    break;
                }
            }
        }

        if (valid_prefix) {
            id_map.try_insert(e.ref_name);

            if (e.mapq >= min_mapq and (not e.is_not_primary())) {
                mappings[e.read_name].emplace(e);
            }
        }
    });
}


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

    BamAlignment e;
    size_t l = 0;

    while (reader.GetNextAlignment(e) ) {
//        cerr << "e.Name" << ' ' << e.Name << '\n';
//        cerr << "e.RefID" << ' ' << e.RefID << '\n';
//        cerr << "int32_t(l)" << ' ' << int32_t(l) << '\n';
//        cerr << "int16_t(e.AlignmentFlag)" << ' ' << int(e.AlignmentFlag) << '\n';
//        cerr << "int8_t(e.MapQuality)" << ' ' << int(e.MapQuality) << '\n';
//        cerr << '\n';

        auto& ref_name = reference_data.at(e.RefID).RefName;

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

            if (e.MapQuality >= min_mapq and e.IsPrimaryAlignment()) {
                SamElement s(e.Name, ref_name, int32_t(l), int16_t(e.AlignmentFlag), int8_t(e.MapQuality));
                mappings[e.Name].emplace(s);
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
        contact_map_t& contact_map
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

                        contact_map[id][id2]++;
                        contact_map[id2][id]++;
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
        const contact_map_t& contact_map,
        const IncrementalIdMap<string>& id_map){
    ofstream output_file(output_path);

    if (not output_file.is_open() or not output_file.good()){
        throw std::runtime_error("ERROR: could not write to file: " + output_path.string());
    }

    for (const auto& [id,map2]: contact_map){
        for (const auto& [id2,count]: map2){
            output_file << id_map.get_name(id) << ',' << id_map.get_name(id2) << ',' << count << '\n';
        }
    }

}


void write_config(
        path output_path,
        path sam_path,
        string required_prefix,
        int8_t min_mapq,
        size_t n_threads){

    ofstream file(output_path);

    if (not file.is_open() or not file.good()){
        throw std::runtime_error("ERROR: could not write to file: " + output_path.string());
    }

    file << "sam_path" << ',' << sam_path << '\n';
    file << "required_prefix" << ',' << required_prefix << '\n';
    file << "min_mapq" << ',' << min_mapq << '\n';
    file << "n_threads" << ',' << n_threads << '\n';
}


void phase_hic(path output_dir, path sam_path, string required_prefix, int8_t min_mapq, size_t n_threads){
    if (exists(output_dir)){
        throw runtime_error("ERROR: output directory exists already");
    }
    else {
        create_directories(output_dir);
    }

    // Id-to-name bimap for reference contigs
    IncrementalIdMap<string> id_map(true);

    // TODO: Load GFA first, and find a method for identifying haplotypic bubbles
    // TODO: in other words, stop relying on shasta conventions

    // Mappings grouped by their read name (to simplify paired-end parsing)
    unpaired_mappings_t mappings;

    // Datastructures to represent linkages from hiC
    contact_map_t contact_map;
    vector <vector <int32_t> > adjacency;

    if (sam_path.extension() == ".sam") {
        parse_unpaired_sam_file(sam_path, mappings, id_map, required_prefix, min_mapq);
    }
    if (sam_path.extension() == ".bam"){
        parse_unpaired_bam_file(sam_path, mappings, id_map, required_prefix, min_mapq);
    }
    else{
        throw runtime_error("ERROR: unrecognized extension for SAM/BAM input file: " + sam_path.extension().string());
    }

    remove_singleton_reads(mappings);

    // Build the contact map by iterating the pairs and creating edges in an all-by-all fashion between pairs
    generate_contact_map_from_mappings(mappings, id_map, contact_map);

    // To keep track of pairs of segments which exist in diploid bubbles
    BubbleGraph bubble_graph(id_map, contact_map);

    cerr << "Phasing " << bubble_graph.size() << " bubbles" << '\n';

    phase_contacts(contact_map, id_map, bubble_graph, n_threads);

    int64_t score = compute_total_consistency_score(bubble_graph, contact_map);

    string suffix1 = "p" + required_prefix;
    string suffix2 = "m" + to_string(int(min_mapq));
    string suffix3 = "s" + to_string(int(score));

    path contacts_output_path = output_dir / "contacts.csv";
    path phases_output_path = output_dir / "phases.csv";
    path config_output_path = output_dir / "config.csv";

    write_contact_map(contacts_output_path, contact_map, id_map);
    bubble_graph.write_bandage_csv(phases_output_path, id_map);
}


void chain_phased_gfa(path gfa_path, IncrementalIdMap<string>& id_map, const BubbleGraph& bubble_graph){
    HashGraph graph;

    gfa_to_handle_graph(graph, id_map, gfa_path, false);

    vector<HashGraph> connected_components;
    vector <IncrementalIdMap<string> > connected_component_ids;

    cerr << "Finding connected components..." << '\n';

    split_connected_components(graph, id_map, connected_components, connected_component_ids, true);

//    for (size_t c=0; c<connected_components.size(); c++){
//        unzip(connected_components[c], connected_component_ids[c], false);
//
//        auto& cc_graph = connected_components[c];
//        auto& cc_id_map = connected_component_ids[c];
//
//        // Generate criteria for diploid node BFS
//        unordered_set<nid_t> diploid_nodes;
//        generate_ploidy_critera(cc_graph, cc_id_map, diploid_names, diploid_nodes);
//
//        Bipartition ploidy_bipartition(cc_graph, cc_id_map, diploid_nodes);
//        ploidy_bipartition.partition();
//
//        // Generate criteria for node-chaining BFS
//        unordered_set<nid_t> chain_nodes;
//        generate_chain_critera(ploidy_bipartition, chain_nodes);
//
//        Bipartition chain_bipartition(cc_graph, cc_id_map, chain_nodes);
//        chain_bipartition.partition();
//
//        string filename_prefix = "component_" + to_string(c) + "_";
//        ofstream file(filename_prefix + ".gfa");
//        handle_graph_to_gfa(connected_components[c], connected_component_ids[c], file);
//
//        ofstream test_gfa_meta(filename_prefix + "ploidy_metagraph.gfa");
//        handle_graph_to_gfa(ploidy_bipartition.metagraph, test_gfa_meta);
//
//        ofstream test_csv_meta_ploidy(filename_prefix + "ploidy_metagraph.csv");
//        ploidy_bipartition.write_meta_graph_csv(test_csv_meta_ploidy);
//
//        ofstream test_csv_parent_ploidy(filename_prefix + "ploidy_parent_graph.csv");
//        ploidy_bipartition.write_parent_graph_csv(test_csv_parent_ploidy);
//
//        ofstream test_gfa_chain(filename_prefix + "chain_metagraph.gfa");
//        handle_graph_to_gfa(chain_bipartition.metagraph, test_gfa_chain);
//
//        ofstream test_csv_meta_chain(filename_prefix + "chain_metagraph.csv");
//        chain_bipartition.write_meta_graph_csv(test_csv_meta_chain);
//
//        ofstream test_csv_parent_chain(filename_prefix + "chain_parent_graph.csv");
//        chain_bipartition.write_parent_graph_csv(test_csv_parent_chain);
//
//        unordered_set<string> paternal_node_names;
//        unordered_set<string> maternal_node_names;
//
//        merge_diploid_singletons(diploid_names, chain_bipartition);
//
//        ofstream test_gfa_chain_merged(filename_prefix + "chain_metagraph_merged.gfa");
//        handle_graph_to_gfa(chain_bipartition.metagraph, test_gfa_chain_merged);
//
//        ofstream test_csv_meta_chain_merged(filename_prefix + "chain_metagraph_merged.csv");
//        chain_bipartition.write_meta_graph_csv(test_csv_meta_chain_merged);
//
//        ofstream test_csv_parent_chain_merged(filename_prefix + "chain_parent_graph_merged.csv");
//        chain_bipartition.write_parent_graph_csv(test_csv_parent_chain_merged);
//
//        string component_path_prefix = to_string(c);
//        phase_chains<T,T2>(
//                chain_bipartition,
//                cc_graph,
//                cc_id_map,
//                diploid_path_names,
//                ks,
//                paternal_node_names,
//                maternal_node_names,
//                provenance_csv_file,
//                component_path_prefix,
//                path_delimiter);
//
//        ofstream test_gfa_phased(filename_prefix + "phased.gfa");
//        handle_graph_to_gfa(cc_graph, cc_id_map, test_gfa_phased);
//
//        cc_graph.for_each_handle([&](const handle_t& h){
//            auto id = cc_graph.get_id(h);
//            auto name = cc_id_map.get_name(id);
//
//            if (maternal_node_names.count(name) > 0){
//                maternal_fasta << '>' << name << '\n';
//                maternal_fasta << cc_graph.get_sequence(h) << '\n';
//            }
//            else if (paternal_node_names.count(name) > 0){
//                paternal_fasta << '>' << name << '\n';
//                paternal_fasta << cc_graph.get_sequence(h) << '\n';
//            }
//            else {
//                unphased_initial_fasta << '>' << name << '\n';
//                unphased_initial_fasta << cc_graph.get_sequence(h) << '\n';
//                unphased_handles_per_component[c].emplace_back(h);
//            }
//        });
//    }

}


int main (int argc, char* argv[]){
    path sam_path;
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

    phase_hic(output_dir, sam_path, required_prefix, min_mapq, n_threads);

    return 0;
}
