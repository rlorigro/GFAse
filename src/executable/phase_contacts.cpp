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
using gfase::Bipartition;
using gfase::BubbleGraph;
using gfase::SamElement;
using gfase::Bubble;

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
        path output_dir,
        path sam_path,
        string required_prefix,
        int8_t min_mapq,
        size_t n_threads){

    path output_path = output_dir / "config.csv";
    ofstream file(output_path);

    if (not file.is_open() or not file.good()){
        throw std::runtime_error("ERROR: could not write to file: " + output_path.string());
    }

    file << "sam_path" << ',' << sam_path << '\n';
    file << "required_prefix" << ',' << required_prefix << '\n';
    file << "min_mapq" << ',' << int(min_mapq) << '\n';
    file << "n_threads" << ',' << n_threads << '\n';
}


void generate_ploidy_criteria_from_bubble_graph(
        const PathHandleGraph& graph,
        const IncrementalIdMap<string>& id_map,
        const BubbleGraph& bubble_graph,
        unordered_set<nid_t>& diploid_nodes
){

    bubble_graph.for_each_node_id([&](const int32_t id){
        // Node names for haplotypes should match the paths that they were created from
        auto name = id_map.get_name(id);

        // Add diploid nodes to the set
        diploid_nodes.emplace(id);

        return true;
    });
}


void merge_diploid_singletons(const BubbleGraph& bubble_graph, Bipartition& chain_bipartition){
    unordered_set <pair <size_t,size_t> > to_be_merged;

    chain_bipartition.for_each_subgraph([&](const HandleGraph& subgraph, size_t subgraph_index, bool partition){

        // Look for phaseable subgraphs only (they contain at least one diploid node) and size == 1 (singleton)
        if (chain_bipartition.get_partition_of_subgraph(subgraph_index) == 0 and chain_bipartition.get_subgraph_size(subgraph_index) == 1){
            nid_t singleton_id;
            string singleton_name;

            chain_bipartition.for_each_handle_in_subgraph(subgraph_index, [&](const handle_t& h){
                singleton_id = chain_bipartition.get_id_of_parent_handle(h);
                singleton_name = chain_bipartition.get_name_of_parent_node(singleton_id);
            });

            // Find other diploid node and verify is also singleton
            nid_t other_id = bubble_graph.get_other_side(int32_t(singleton_id));

            auto other_subgraph_index = chain_bipartition.get_subgraph_index_of_parent_node(other_id);

            // Use a defined ordering of singleton pairs to keep track of which have been visited
            to_be_merged.emplace(min(subgraph_index,other_subgraph_index), max(subgraph_index,other_subgraph_index));
        }
    });

    for (auto& item: to_be_merged){
        chain_bipartition.merge_subgraphs(item.first, item.second);
    }
}


void write_chaining_info_to_file(
        path output_dir,
        const Bipartition& ploidy_bipartition,
        const Bipartition& chain_bipartition,
        const PathHandleGraph& graph,
        const IncrementalIdMap<string>& id_map,
        const string& filename_prefix,
        size_t component_index
        ){

    size_t c = component_index;

    path file_path = output_dir / "components" / to_string(c) / (filename_prefix + ".gfa");
    ofstream file(file_path);
    handle_graph_to_gfa(graph, id_map, file);

    path test_gfa_chain_path = output_dir / "components" / to_string(c) / (filename_prefix + "chain_metagraph.gfa");
    ofstream test_gfa_chain(test_gfa_chain_path);

    if (not test_gfa_chain.is_open() or not test_gfa_chain.good()){
        throw runtime_error("ERROR: could not write to file: " + test_gfa_chain_path.string());
    }
    handle_graph_to_gfa(chain_bipartition.metagraph, test_gfa_chain);

    path test_csv_meta_chain_path = output_dir / "components" / to_string(c) / (filename_prefix + "chain_metagraph.csv");
    ofstream test_csv_meta_chain(test_csv_meta_chain_path);

    if (not test_csv_meta_chain.is_open() or not test_csv_meta_chain.good()){
        throw runtime_error("ERROR: could not write to file: " + test_csv_meta_chain_path.string());
    }
    chain_bipartition.write_meta_graph_csv(test_csv_meta_chain);

    path test_csv_parent_chain_path = output_dir / "components" / to_string(c) / (filename_prefix + "chain_parent_graph.csv");
    ofstream test_csv_parent_chain(test_csv_parent_chain_path);

    if (not test_csv_parent_chain.is_open() or not test_csv_parent_chain.good()){
        throw runtime_error("ERROR: could not write to file: " + test_csv_parent_chain_path.string());
    }
    chain_bipartition.write_parent_graph_csv(test_csv_parent_chain);

    path test_gfa_chain_merged_path = output_dir / "components" / to_string(c) / (filename_prefix + "chain_metagraph_merged.gfa");
    ofstream test_gfa_chain_merged(test_gfa_chain_merged_path);

    if (not test_gfa_chain_merged.is_open() or not test_gfa_chain_merged.good()){
        throw runtime_error("ERROR: could not write to file: " + test_gfa_chain_merged_path.string());
    }
    handle_graph_to_gfa(chain_bipartition.metagraph, test_gfa_chain_merged);

    path test_csv_meta_chain_merged_path = output_dir / "components" / to_string(c) / (filename_prefix + "chain_metagraph_merged.csv");
    ofstream test_csv_meta_chain_merged(test_csv_meta_chain_merged_path);

    if (not test_csv_meta_chain_merged.is_open() or not test_csv_meta_chain_merged.good()){
        throw runtime_error("ERROR: could not write to file: " + test_csv_meta_chain_merged_path.string());
    }
    chain_bipartition.write_meta_graph_csv(test_csv_meta_chain_merged);

    path test_csv_parent_chain_merged_path = output_dir / "components" / to_string(c) / (filename_prefix + "chain_parent_graph_merged.csv");
    ofstream test_csv_parent_chain_merged(test_csv_parent_chain_merged_path);

    if (not test_csv_parent_chain_merged.is_open() or not test_csv_parent_chain_merged.good()){
        throw runtime_error("ERROR: could not write to file: " + test_csv_parent_chain_merged_path.string());
    }
    chain_bipartition.write_parent_graph_csv(test_csv_parent_chain_merged);
}


void chain_phased_gfa(path gfa_path, IncrementalIdMap<string>& id_map, const BubbleGraph& bubble_graph, path output_dir){
    HashGraph graph;

    // When assigning IDs to new nodes in the graph, gfa_to_handle should reuse existing ones, so there
    // will be no conflicts in the BubbleGraph IDs and the graph IDs (as long as id_map is 1-based and size_t)
    gfa_to_handle_graph(graph, id_map, gfa_path, false);

    vector<HashGraph> connected_components;
    vector <IncrementalIdMap<string> > connected_component_ids;

    cerr << "Finding connected components..." << '\n';

    split_connected_components(graph, id_map, connected_components, true);

    cerr << "Chaining phased bubbles..." << '\n';

    for (size_t c=0; c<connected_components.size(); c++){
        unzip(connected_components[c], connected_component_ids[c], false);

        auto& cc_graph = connected_components[c];

        // Generate criteria for diploid node BFS
        unordered_set<nid_t> diploid_nodes;
        generate_ploidy_criteria_from_bubble_graph(cc_graph, id_map, bubble_graph, diploid_nodes);

        Bipartition ploidy_bipartition(cc_graph, id_map, diploid_nodes);
        ploidy_bipartition.partition();

        // Generate criteria for node-chaining BFS
        unordered_set<nid_t> chain_nodes;
        generate_chain_critera(ploidy_bipartition, chain_nodes);

        Bipartition chain_bipartition(cc_graph, id_map, chain_nodes);
        chain_bipartition.partition();

        create_directories(output_dir / "components" / to_string(c));

        unordered_set<string> phase_0_node_names;
        unordered_set<string> phase_1_node_names;

        merge_diploid_singletons(bubble_graph, chain_bipartition);

        string filename_prefix = "component_" + to_string(c) + "_";

        path test_gfa_meta_path = output_dir / "components" / to_string(c) / (filename_prefix + "ploidy_metagraph.gfa");
        ofstream test_gfa_meta(test_gfa_meta_path);

        if (not test_gfa_meta.is_open() or not test_gfa_meta.good()){
            throw runtime_error("ERROR: could not write to file: " + test_gfa_meta_path.string());
        }
        handle_graph_to_gfa(ploidy_bipartition.metagraph, test_gfa_meta);

        path test_csv_meta_ploidy_path = output_dir / "components" / to_string(c) / (filename_prefix + "ploidy_metagraph.csv");
        ofstream test_csv_meta_ploidy(test_csv_meta_ploidy_path);

        if (not test_csv_meta_ploidy.is_open() or not test_csv_meta_ploidy.good()){
            throw runtime_error("ERROR: could not write to file: " + test_csv_meta_ploidy_path.string());
        }
        ploidy_bipartition.write_meta_graph_csv(test_csv_meta_ploidy);

        path test_csv_parent_ploidy_path = output_dir / "components" / to_string(c) / (filename_prefix + "ploidy_parent_graph.csv");
        ofstream test_csv_parent_ploidy(test_csv_parent_ploidy_path);

        if (not test_csv_parent_ploidy.is_open() or not test_csv_parent_ploidy.good()){
            throw runtime_error("ERROR: could not write to file: " + test_csv_parent_ploidy_path.string());
        }
        ploidy_bipartition.write_parent_graph_csv(test_csv_parent_ploidy);

        write_chaining_info_to_file(output_dir, ploidy_bipartition, chain_bipartition, cc_graph, id_map, filename_prefix, c);

        path_handle_t phase_0_path;
        path_handle_t phase_1_path;

        path phase_0_fasta_path = output_dir / "phase_0.fasta";
        ofstream phase_0_fasta(phase_0_fasta_path);
        path phase_1_fasta_path = output_dir / "phase_1.fasta";
        ofstream phase_1_fasta(phase_1_fasta_path);
        path unphased_initial_fasta_path = output_dir / "unphased_initial.fasta";
        ofstream unphased_initial_fasta(unphased_initial_fasta_path);
        path unphased_fasta_path = output_dir / "unphased.fasta";
        ofstream unphased_fasta(unphased_fasta_path);
        path provenance_csv_file_path = output_dir / "phase_chains.csv";
        ofstream provenance_csv_file(provenance_csv_file_path);

        if (not (phase_0_fasta.is_open() and provenance_csv_file.good())){
            throw runtime_error("ERROR: file could not be written: " + phase_0_fasta_path.string());
        }

        if (not (phase_1_fasta.is_open() and provenance_csv_file.good())){
            throw runtime_error("ERROR: file could not be written: " + phase_1_fasta_path.string());
        }

        if (not (unphased_initial_fasta.is_open() and provenance_csv_file.good())){
            throw runtime_error("ERROR: file could not be written: " + unphased_initial_fasta_path.string());
        }

        if (not (unphased_fasta.is_open() and provenance_csv_file.good())){
            throw runtime_error("ERROR: file could not be written: " + unphased_fasta_path.string());
        }

        if (not (provenance_csv_file.is_open() and provenance_csv_file.good())){
            throw runtime_error("ERROR: file could not be written: " + provenance_csv_file_path.string());
        }

        provenance_csv_file << "path_name" << ',' << "n_steps" << ',' << "nodes" << '\n';

        vector <vector <handle_t> > unphased_handles_per_component(connected_components.size());

        string component_path_prefix = to_string(c);
        auto prev_subgraph_index = numeric_limits<size_t>::max();
        for_element_in_bubble_chain(
                chain_bipartition,
                cc_graph,
                id_map,
                bubble_graph,
                [&](const vector<string>& node_names, size_t subgraph_index){

                    if (subgraph_index != prev_subgraph_index){
                        string path_prefix = component_path_prefix + '.' + to_string(subgraph_index);

                        string phase_0_path_name = path_prefix + ".0";
                        string phase_1_path_name = path_prefix + ".1";

                        phase_0_path = cc_graph.create_path_handle(phase_0_path_name);
                        phase_1_path = cc_graph.create_path_handle(phase_1_path_name);

                        phase_0_node_names.emplace(phase_0_path_name);
                        phase_1_node_names.emplace(phase_1_path_name);
                    }

                    // If there's only one node in the element, then it must be in-between bubbles
                    if (node_names.size() == 1){
                        auto node = cc_graph.get_handle(id_map.get_id(node_names[0]));

                        cc_graph.append_step(phase_0_path, node);
                        cc_graph.append_step(phase_1_path, node);
                    }
                    // If there are 2 nodes then it must be a bubble
                    else{
                        auto& name_a = node_names[0];
                        auto& name_b = node_names[1];

                        auto id_a = id_map.get_id(name_a);
                        auto id_b = id_map.get_id(name_b);

                        auto b = bubble_graph.get_bubble_of_node(int32_t(id_a));

                        auto node_a = cc_graph.get_handle(id_a);
                        auto node_b = cc_graph.get_handle(id_b);

                        // Choose the more supported orientation, defaulting to "forward orientation" if equal
                        if (b.phase == 0){
                            cc_graph.append_step(phase_0_path, node_a);
                            cc_graph.append_step(phase_1_path, node_b);
                        }
                        else{
                            cc_graph.append_step(phase_1_path, node_a);
                            cc_graph.append_step(phase_0_path, node_b);
                        }
                    }

                    prev_subgraph_index = subgraph_index;
                });

        unzip(cc_graph, id_map, false);
        write_paths_to_csv(cc_graph, id_map, provenance_csv_file);

        path test_gfa_phased_path = output_dir / "components" / to_string(c) / (filename_prefix + "phased.gfa");
        ofstream test_gfa_phased(test_gfa_phased_path);
        if (not test_gfa_phased.is_open() or not test_gfa_phased.good()){
            throw runtime_error("ERROR: could not write to file: " + test_gfa_phased_path.string());
        }

        handle_graph_to_gfa(cc_graph, id_map, test_gfa_phased);

        cc_graph.for_each_handle([&](const handle_t& h){
            auto id = cc_graph.get_id(h);
            auto name = id_map.get_name(id);

            cerr << c << ' ' << id << ' ' << name;

            if (phase_0_node_names.count(name) > 0){
                cerr << ' ' << '0' << '\n';
                phase_0_fasta << '>' << name << '\n';
                phase_0_fasta << cc_graph.get_sequence(h) << '\n';
            }
            else if (phase_1_node_names.count(name) > 0){
                cerr << ' ' << '1' << '\n';
                phase_1_fasta << '>' << name << '\n';
                phase_1_fasta << cc_graph.get_sequence(h) << '\n';
            }
            else {
                cerr << ' ' << "unphased" << '\n';
                unphased_initial_fasta << '>' << name << '\n';
                unphased_initial_fasta << cc_graph.get_sequence(h) << '\n';
                unphased_handles_per_component[c].emplace_back(h);
            }
        });
    }
}


void phase_hic(path output_dir, path sam_path, path gfa_path, string required_prefix, int8_t min_mapq, size_t n_threads){
    if (exists(output_dir)){
        throw runtime_error("ERROR: output directory exists already");
    }
    else {
        create_directories(output_dir);
    }

    write_config(output_dir, sam_path, required_prefix, min_mapq, n_threads);

    // Id-to-name bimap for reference contigs
    IncrementalIdMap<string> id_map(false);

    // TODO: Load GFA first, and find a method for identifying haplotypic bubbles
    // TODO: in other words, stop relying on shasta conventions

    // Mappings grouped by their read name (to simplify paired-end parsing)
    unpaired_mappings_t mappings;

    // Datastructures to represent linkages from hiC
    contact_map_t contact_map;
    vector <vector <int32_t> > adjacency;

    cerr << "Loading alignments..." << '\n';

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

    cerr << "Generating contact map..." << '\n';

    // Build the contact map by iterating sets of alignments and creating edges in an all-by-all fashion within sets
    generate_contact_map_from_mappings(mappings, id_map, contact_map);

    cerr << "Constructing bubble graph..." << '\n';

    // TODO: replace name-based bubble finding with sketch or alignment-based
    // To keep track of pairs of segments which exist in diploid bubbles
    BubbleGraph bubble_graph(id_map, contact_map);

    cerr << "Phasing " << bubble_graph.size() << " bubbles" << '\n';

    phase_contacts(contact_map, id_map, bubble_graph, n_threads);

    int64_t score = compute_total_consistency_score(bubble_graph, contact_map);

    string suffix1 = "p" + required_prefix;
    string suffix2 = "m" + to_string(int(min_mapq));
    string suffix3 = "s" + to_string(int(score));

    cerr << "Writing phasing results to disk... " << '\n';

    path contacts_output_path = output_dir / "contacts.csv";
    path phases_output_path = output_dir / "phases.csv";
    path config_output_path = output_dir / "config.csv";

    write_contact_map(contacts_output_path, contact_map, id_map);
    bubble_graph.write_bandage_csv(phases_output_path, id_map);

    if (not gfa_path.empty()){
        chain_phased_gfa(gfa_path, id_map, bubble_graph, output_dir);
    }
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

    phase_hic(output_dir, sam_path, gfa_path, required_prefix, min_mapq, n_threads);

    return 0;
}
