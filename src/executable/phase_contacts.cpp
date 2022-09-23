#include "IncrementalIdMap.hpp"
#include "gfa_to_handle.hpp"
#include "graph_utility.hpp"
#include "BubbleGraph.hpp"
#include "Bipartition.hpp"
#include "hash_graph.hpp"
#include "Filesystem.hpp"
#include "Timer.hpp"
#include "CLI11.hpp"
#include "Sam.hpp"
#include "Bam.hpp"

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
using gfase::Timer;
using gfase::Bam;

using bdsg::HashGraph;


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


using weighted_contact_map_t = sparse_hash_map <int32_t, sparse_hash_map<int32_t, map <uint8_t, int32_t> > >;


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


void update_contact_map(
        vector<SamElement>& alignments,
        contact_map_t& contact_map,
        IncrementalIdMap<string>& id_map){

    // Iterate one triangle of the all-by-all matrix, adding up mapqs for reads on both end of the pair
    for (size_t i=0; i<alignments.size(); i++){
//        cerr << alignments[i] << '\n';

        for (size_t j=i+1; j<alignments.size(); j++) {
            auto& a = alignments[i];
            auto& b = alignments[j];

            auto ref_id_a = id_map.try_insert(a.ref_name);
            auto ref_id_b = id_map.try_insert(b.ref_name);

            // TODO: split left and right mapq instead of taking min?
            contact_map[ref_id_a][ref_id_b]++;
            contact_map[ref_id_b][ref_id_a]++;
        }
    }
}


void parse_unpaired_bam_file(
        path bam_path,
        contact_map_t& contact_map,
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


void chain_phased_gfa(MutablePathDeletableHandleGraph& graph, IncrementalIdMap<string>& id_map, const BubbleGraph& bubble_graph, path output_dir){

    vector<HashGraph> connected_components;
    vector <IncrementalIdMap<string> > connected_component_ids;

    cerr << "Finding connected components..." << '\n';

    split_connected_components(graph, id_map, connected_components, true);

    cerr << "Chaining phased bubbles..." << '\n';

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

    if (not (phase_0_fasta.is_open() and phase_0_fasta.good())){
        throw runtime_error("ERROR: file could not be written: " + phase_0_fasta_path.string());
    }

    if (not (phase_1_fasta.is_open() and phase_1_fasta.good())){
        throw runtime_error("ERROR: file could not be written: " + phase_1_fasta_path.string());
    }

    if (not (unphased_initial_fasta.is_open() and unphased_initial_fasta.good())){
        throw runtime_error("ERROR: file could not be written: " + unphased_initial_fasta_path.string());
    }

    if (not (unphased_fasta.is_open() and unphased_fasta.good())){
        throw runtime_error("ERROR: file could not be written: " + unphased_fasta_path.string());
    }

    if (not (provenance_csv_file.is_open() and provenance_csv_file.good())){
        throw runtime_error("ERROR: file could not be written: " + provenance_csv_file_path.string());
    }

    provenance_csv_file << "path_name" << ',' << "n_steps" << ',' << "nodes" << '\n';

    for (size_t c=0; c<connected_components.size(); c++){
        cerr << c << '\n';
        unzip(connected_components[c], connected_component_ids[c], false);

        auto& cc_graph = connected_components[c];

        cerr << "Generating diploid criteria" << '\n';

        // Generate criteria for diploid node BFS
        unordered_set<nid_t> diploid_nodes;
        generate_ploidy_criteria_from_bubble_graph(cc_graph, id_map, bubble_graph, diploid_nodes);

        cerr << "Doing ploidy partition" << '\n';

        Bipartition ploidy_bipartition(cc_graph, id_map, diploid_nodes);
        ploidy_bipartition.partition();

        cerr << "Generating chaining criteria" << '\n';

        // Generate criteria for node-chaining BFS
        unordered_set<nid_t> chain_nodes;
        generate_chain_critera(ploidy_bipartition, chain_nodes);

        cerr << "Doing chain partition" << '\n';

        Bipartition chain_bipartition(cc_graph, id_map, chain_nodes);
        chain_bipartition.partition();

        create_directories(output_dir / "components" / to_string(c));

        unordered_set<string> phase_0_node_names;
        unordered_set<string> phase_1_node_names;

        cerr << "Merging singletons" << '\n';

        merge_diploid_singletons(bubble_graph, chain_bipartition);

        cerr << "Writing results" << '\n';

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

                        if (bubble_graph.find_bubble_id_of_node(int32_t(id_a)) != bubble_graph.find_bubble_id_of_node(int32_t(id_b))){
                            throw runtime_error("ERROR: nodes in diploid portion of chain are not of same bubble in bubble graph: " + (name_a + ',' + name_b));
                        }

                        // Only need one id to reach the corresponding bubble
                        auto bubble = bubble_graph.get_bubble_of_node(int32_t(id_a));

                        nid_t id_0 = bubble.first();
                        nid_t id_1 = bubble.second();

                        auto node_0 = cc_graph.get_handle(id_0);
                        auto node_1 = cc_graph.get_handle(id_1);

                        // Choose orientation supported by hi-c linkage (decided already, during phasing)
                        cc_graph.append_step(phase_0_path, node_0);
                        cc_graph.append_step(phase_1_path, node_1);
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

            if (phase_0_node_names.count(name) > 0){
                phase_0_fasta << '>' << name << '\n';
                phase_0_fasta << cc_graph.get_sequence(h) << '\n';
            }
            else if (phase_1_node_names.count(name) > 0){
                phase_1_fasta << '>' << name << '\n';
                phase_1_fasta << cc_graph.get_sequence(h) << '\n';
            }
            else {
                unphased_initial_fasta << '>' << name << '\n';
                unphased_initial_fasta << cc_graph.get_sequence(h) << '\n';
                unphased_handles_per_component[c].emplace_back(h);
            }
        });
    }
}


void phase_hic(path output_dir, path sam_path, path gfa_path, string required_prefix, int8_t min_mapq, size_t n_threads){
    Timer t;

    if (exists(output_dir)){
        throw runtime_error("ERROR: output directory exists already");
    }
    else {
        create_directories(output_dir);
    }

    write_config(output_dir, sam_path, required_prefix, min_mapq, n_threads);

    // Id-to-name bimap for reference contigs
    IncrementalIdMap<string> id_map(false);

    // Datastructures to represent linkages from hiC
    contact_map_t contact_map;
    vector <vector <int32_t> > adjacency;

    cerr << t << "Loading alignments as contact map..." << '\n';

    if (sam_path.extension() == ".bam"){
        parse_unpaired_bam_file(sam_path, contact_map, id_map, required_prefix, min_mapq);
    }
    else{
        throw runtime_error("ERROR: unrecognized extension for BAM input file: " + sam_path.extension().string());
    }

    HashGraph graph;

    // To keep track of pairs of segments which exist in diploid bubbles
    BubbleGraph bubble_graph;

    if (gfa_path.empty()) {
        cerr << t << "Constructing bubble graph..." << '\n';

        // Initialize using shasta naming convention
        bubble_graph = BubbleGraph(id_map, contact_map);
    }
    else{
        cerr << t << "GFA provided - Loading graph..." << '\n';

        // Construct graph from GFA
        gfa_to_handle_graph(graph, id_map, gfa_path, false);

        cerr << t << "Constructing bubble graph..." << '\n';

        // Initialize using shasta naming convention
        bubble_graph = BubbleGraph(id_map, contact_map);
    }

    cerr << t << "Phasing " << bubble_graph.size() << " bubbles" << '\n';

    phase_contacts(contact_map, id_map, bubble_graph, n_threads);

    int64_t score = compute_total_consistency_score(bubble_graph, contact_map);

    string suffix1 = "p" + required_prefix;
    string suffix2 = "m" + to_string(int(min_mapq));
    string suffix3 = "s" + to_string(int(score));

    cerr << t << "Writing phasing results to file... " << '\n';

    path contacts_output_path = output_dir / "contacts.csv";
    path phases_output_path = output_dir / "phases.csv";
    path config_output_path = output_dir / "config.csv";

    write_contact_map(contacts_output_path, contact_map, id_map);
    bubble_graph.write_bandage_csv(phases_output_path, id_map);

    if (not gfa_path.empty()){
        chain_phased_gfa(graph, id_map, bubble_graph, output_dir);
    }

    cerr << t << "Done" << '\n';
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
