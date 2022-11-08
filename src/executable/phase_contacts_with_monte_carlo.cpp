#include "MultiContactGraph.hpp"
#include "IncrementalIdMap.hpp"
#include "gfa_to_handle.hpp"
#include "handle_to_gfa.hpp"
#include "graph_utility.hpp"
#include "hash_graph.hpp"
#include "Filesystem.hpp"
#include "optimize.hpp"
#include "Chainer.hpp"
#include "Timer.hpp"
#include "CLI11.hpp"
#include "Sam.hpp"
#include "Bam.hpp"

using gfase::for_element_in_sam_file;
using gfase::unpaired_mappings_t;
using gfase::paired_mappings_t;
using gfase::contact_map_t;
using gfase::unzip;

using gfase::gfa_to_handle_graph;
using gfase::handle_graph_to_gfa;
using gfase::MultiContactGraph;
using gfase::IncrementalIdMap;
using gfase::SamElement;
using gfase::Chainer;
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
        MultiContactGraph& contact_graph,
        IncrementalIdMap<string>& id_map,
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

        // Only allow reads with mapq > min_mapq and not secondary
        if (a.mapq >= min_mapq and a.is_primary()) {
            alignments.emplace_back(a);
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
        int8_t min_mapq,
        size_t n_threads){

    path output_path = output_dir / "config.csv";
    ofstream file(output_path);

    if (not file.is_open() or not file.good()){
        throw std::runtime_error("ERROR: could not write to file: " + output_path.string());
    }

    file << "sam_path" << ',' << sam_path << '\n';
    file << "min_mapq" << ',' << int(min_mapq) << '\n';
    file << "n_threads" << ',' << n_threads << '\n';
}


void write_gfa_to_file(PathHandleGraph& graph, IncrementalIdMap<string>& id_map, path output_gfa_path){
    ofstream chained_gfa(output_gfa_path);

    if (not (chained_gfa.is_open() and chained_gfa.good())){
        throw runtime_error("ERROR: could not write to file: " + output_gfa_path.string());
    }

    handle_graph_to_gfa(graph, id_map, chained_gfa);

}


void write_nodes_to_fasta(
        const HandleGraph& graph,
        const IncrementalIdMap<string>& id_map,
        const MultiContactGraph& contact_graph,
        const Chainer& chainer,
        path output_dir
        ){

    path phase_0_fasta_path = output_dir / "phase_0.fasta";
    path phase_1_fasta_path = output_dir / "phase_1.fasta";
    path unphased_fasta_path = output_dir / "unphased.fasta";

    ofstream phase_0_fasta(phase_0_fasta_path);
    ofstream phase_1_fasta(phase_1_fasta_path);
    ofstream unphased_fasta(unphased_fasta_path);

    if (not (phase_0_fasta.is_open() and phase_0_fasta.good())){
        throw runtime_error("ERROR: file could not be written: " + phase_0_fasta_path.string());
    }

    if (not (phase_1_fasta.is_open() and phase_1_fasta.good())){
        throw runtime_error("ERROR: file could not be written: " + phase_1_fasta_path.string());
    }

    if (not (unphased_fasta.is_open() and unphased_fasta.good())){
        throw runtime_error("ERROR: file could not be written: " + unphased_fasta_path.string());
    }

    graph.for_each_handle([&](const handle_t& h){
        auto id = graph.get_id(h);
        auto name = id_map.get_name(id);

        bool in_contact_graph = contact_graph.has_node(int32_t(id));
        bool in_phase_chains = chainer.has_phase_chain(name);

        int8_t partition;

        // Nodes may have been deleted during unzipping, so check the chainer and the contact_graph for their phase
        if (in_contact_graph and not in_phase_chains){
            partition = contact_graph.get_partition(int32_t(id));
        }
        else if (in_phase_chains and not in_contact_graph){
            partition = chainer.get_partition(name);
        }
        else if (not in_contact_graph and not in_phase_chains){
            partition = 0;
        }
        else{
            throw runtime_error("ERROR: node in both phase chains and contact graph: " + name);
        }

        if (partition == -1){
            phase_0_fasta << '>' << name << '\n';
            phase_0_fasta << graph.get_sequence(h) << '\n';
        }
        else if (partition == 1){
            phase_1_fasta << '>' << name << '\n';
            phase_1_fasta << graph.get_sequence(h) << '\n';
        }
        else{
            unphased_fasta << '>' << name << '\n';
            unphased_fasta << graph.get_sequence(h) << '\n';
        }
    });
}


void phase_hic(path output_dir, path sam_path, path gfa_path, int8_t min_mapq, size_t n_threads){
    Timer t;

    if (exists(output_dir)){
        throw runtime_error("ERROR: output directory exists already");
    }
    else {
        create_directories(output_dir);
    }

    path config_output_path = output_dir / "config.csv";
    path id_csv_path = output_dir / "ids.csv";
    path contacts_output_path = output_dir / "contacts.csv";
    path phases_output_path = output_dir / "phases.csv";
    path chained_gfa_path = output_dir / "chained.gfa";
    path unzipped_gfa_path = output_dir / "unzipped.gfa";

    write_config(output_dir, sam_path, min_mapq, n_threads);

    // Id-to-name bimap for reference contigs
    IncrementalIdMap<string> id_map(false);

    // How GFA is stored in memory
    HashGraph graph;

    // To keep track of pairs of segments which exist in diploid bubbles
    MultiContactGraph contact_graph;

    // For finding and unzipping bubble chains
    Chainer chainer;

    cerr << t << "Loading GFA..." << '\n';

    // Construct graph from GFA
    gfa_to_handle_graph(graph, id_map, gfa_path, false, true);

    cerr << t << "Writing IDs to file..." << '\n';

    id_map.write_to_csv(id_csv_path);

    cerr << t << "Loading alignments as contact map..." << '\n';

    if (sam_path.extension() == ".bam"){
        parse_unpaired_bam_file(sam_path, contact_graph, id_map, min_mapq);
    }
    else{
        throw runtime_error("ERROR: unrecognized extension for BAM input file: " + sam_path.extension().string());
    }

    cerr << t << "Processing contacts..." << '\n';

    contact_graph.get_alts_from_shasta_names(id_map);

    // Remove nodes that don't have any involvement in bubbles
    vector<int32_t> to_be_deleted;
    contact_graph.for_each_node([&](int32_t id){
        if (not contact_graph.has_alt(id)){
            to_be_deleted.emplace_back(id);
        }
    });

    for (auto& id: to_be_deleted){
        contact_graph.remove_node(id);
    }

    cerr << t << "Optimizing phases..." << '\n';

    contact_graph.write_contact_map(contacts_output_path, id_map);

    // Remove self edges in contact graph (now that they have been written to disk)
    contact_graph.for_each_node([&](int32_t id){
        contact_graph.remove_edge(id,id);
    });

    monte_carlo_phase_contacts(contact_graph, id_map, output_dir, n_threads);

    cerr << t << "Writing phasing results to file... " << '\n';

    contact_graph.write_contact_map(contacts_output_path, id_map);
    contact_graph.write_bandage_csv(phases_output_path, id_map);

    chainer.generate_chain_paths(graph, id_map, contact_graph);
    chainer.write_chainable_nodes_to_bandage_csv(output_dir, id_map);

    write_gfa_to_file(graph, id_map, chained_gfa_path);

    unzip(graph, id_map, false, false);

    cerr << t << "Writing GFA... " << '\n';

    write_gfa_to_file(graph, id_map, unzipped_gfa_path);

    cerr << t << "Writing FASTA... " << '\n';

    write_nodes_to_fasta(graph, id_map, contact_graph, chainer, output_dir);

    cerr << t << "Done" << '\n';
}


int main (int argc, char* argv[]){
    path sam_path;
    path gfa_path;
    path output_dir;
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
            "-m,--min_mapq",
            min_mapq,
            "Minimum required mapq value for mapping to be counted");

    app.add_option(
            "-t,--threads",
            n_threads,
            "Maximum number of threads to use");

    CLI11_PARSE(app, argc, argv);

    phase_hic(output_dir, sam_path, gfa_path, min_mapq, n_threads);

    return 0;
}
