#include "IncrementalIdMap.hpp"
#include "gfa_to_handle.hpp"
#include "graph_utility.hpp"
#include "BubbleGraph.hpp"
#include "Bipartition.hpp"
#include "hash_graph.hpp"
#include "Filesystem.hpp"
#include "Timer.hpp"
#include "CLI11.hpp"
#include "chain.hpp"
#include "Sam.hpp"
#include "Bam.hpp"

using gfase::for_element_in_sam_file;
using gfase::unpaired_mappings_t;
using gfase::paired_mappings_t;
using gfase::contact_map_t;

using gfase::gfa_to_handle_graph;
using gfase::IncrementalIdMap;
using gfase::chain_phased_gfa;
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


void chain(path output_dir, path gfa_path, size_t n_threads, bool use_topology, bool write_sequence){
    Timer t;

    if (exists(output_dir)){
        throw runtime_error("ERROR: output directory exists already");
    }
    else {
        create_directories(output_dir);
    }

    // Id-to-name bimap for reference contigs
    IncrementalIdMap<string> id_map(false);

    HashGraph graph;
    Overlaps overlaps;

    // To keep track of pairs of segments which exist in diploid bubbles
    BubbleGraph bubble_graph;

    if (gfa_path.empty()) {
        cerr << t << "Constructing bubble graph..." << '\n';

        // Initialize using shasta naming convention
        bubble_graph = BubbleGraph(id_map);
    }
    else{
        cerr << t << "GFA provided - Loading graph..." << '\n';

        // Construct graph from GFA
        gfa_to_handle_graph(graph, id_map, overlaps, gfa_path, false);

        cerr << t << "Constructing bubble graph..." << '\n';

        if (use_topology){
            // Initialize using shasta basic topology
            bubble_graph = BubbleGraph(graph);
        }
        else{
            // Initialize using shasta naming convention
            bubble_graph = BubbleGraph(id_map);
        }
    }

    cerr << t << "Writing phasing results to file... " << '\n';

    path phases_output_path = output_dir / "phases.csv";
    path config_output_path = output_dir / "config.csv";

    bubble_graph.write_bandage_csv(phases_output_path, id_map);

    if (not gfa_path.empty()){
        chain_phased_gfa(graph, id_map, overlaps, bubble_graph, output_dir, write_sequence, write_sequence);
    }

    cerr << t << "Done" << '\n';
}


int main (int argc, char* argv[]){
    path gfa_path;
    path output_dir;
    size_t n_threads = 1;
    bool use_topology = false;
    bool write_sequence = false;

    CLI::App app{"App description"};

    app.add_option(
            "-g,--gfa",
            gfa_path,
            "Path to GFA containing assembly graph to be phased")
            ->required();

    app.add_option(
            "-o,--output_dir",
            output_dir,
            "Path to (nonexistent) directory where output will be stored")
            ->required();

    app.add_option(
            "-t,--threads",
            n_threads,
            "Maximum number of threads to use");

    app.add_flag(
        "--use_topology",
        use_topology,
        "Use topology instead of shasta naming convention");

    app.add_flag(
        "--write_sequence",
        write_sequence,
        "Invoke this argument if chained Fasta and GFA are desired");

    CLI11_PARSE(app, argc, argv);

    chain(output_dir, gfa_path, n_threads, use_topology, write_sequence);

    return 0;
}
