#include "IncrementalIdMap.hpp"
#include "gfa_to_handle.hpp"
#include "graph_utility.hpp"
#include "ContactGraph.hpp"
#include "Bipartition.hpp"
#include "hash_graph.hpp"
#include "Filesystem.hpp"
#include "Sequence.hpp"
#include "Hasher2.hpp"
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
using gfase::ContactGraph;
using gfase::Node;
using gfase::Bipartition;
using gfase::SamElement;
using gfase::Sequence;
using gfase::Hasher2;
using gfase::Bubble;
using gfase::Timer;
using gfase::Bam;

using bdsg::HashGraph;
using ghc::filesystem::path;
using CLI::App;


void write_contact_map(
        path output_path,
        const ContactGraph& contact_graph,
        const IncrementalIdMap<string>& id_map){
    ofstream output_file(output_path);

    if (not output_file.is_open() or not output_file.good()){
        throw std::runtime_error("ERROR: could not write to file: " + output_path.string());
    }

    contact_graph.for_each_edge([&](const pair<int32_t,int32_t> edge, int32_t weight){
        output_file << id_map.get_name(edge.first) << ',' << id_map.get_name(edge.second) << ',' << weight << '\n';
    });
}


void write_graph_data(
        path output_directory,
        const Hasher2& hasher,
        const ContactGraph& contact_graph,
        const IncrementalIdMap<string>& id_map){

    path edges_path = output_directory / "edges.csv";
    path nodes_path = output_directory / "nodes.csv";

    ofstream edges_file(edges_path);
    ofstream nodes_file(nodes_path);

    if (not edges_file.is_open() or not edges_file.good()){
        throw std::runtime_error("ERROR: could not write to file: " + edges_path.string());
    }

    if (not nodes_file.is_open() or not nodes_file.good()){
        throw std::runtime_error("ERROR: could not write to file: " + nodes_path.string());
    }

    ContactGraph hash_graph;

    hasher.convert_to_contact_graph(hash_graph, id_map, 0.1, 10, 10);

    // TODO: Renumber nodes from 0-n

    edges_file << "id_a" << ',' << "id_b" << ',' << "name_a" << ',' << "name_b" << ',' << "contact_weight" << ',' << "hash_weight" << '\n';

    // Iterate all contact graph edges and cross-reference with hash graph edges
    contact_graph.for_each_edge([&](const pair<int32_t,int32_t> edge, int32_t weight){
        auto& id_a = edge.first;
        auto& id_b = edge.second;

        // Edges will be ordered by their IDs {min(a,b), max(a,b)}
        auto name_a = id_map.get_name(id_a);
        auto name_b = id_map.get_name(id_b);

        auto hash_weight = hash_graph.get_edge_weight(id_a, id_b);
        auto contact_weight = contact_graph.get_edge_weight(id_a, id_b);

        edges_file << id_a << ',' << id_b << ',' << name_a << ',' << name_b << ',' << contact_weight << ',' << hash_weight << '\n';
    });

    nodes_file << "id" << ',' << "name" << ',' << "length" << ',' << "contact_coverage" << ',' << "hash_coverage" << '\n';

    contact_graph.for_each_node([&](int32_t id, const Node& n){
        auto name = id_map.get_name(id);

        int64_t hash_coverage = 0;
        if (hash_graph.has_node(id)){
            hash_coverage = hash_graph.get_node_coverage(id);
        }

        nodes_file << id << ',' << name << ',' << n.length << ',' << n.coverage << ',' << hash_coverage << '\n';
    });
}


void update_contact_map(
        vector<SamElement>& alignments,
        ContactGraph& contact_graph,
        IncrementalIdMap<string>& id_map){

    // Iterate one triangle of the all-by-all matrix, adding up mapqs for reads on both end of the pair
    for (size_t i=0; i<alignments.size(); i++){
        auto& a = alignments[i];
        auto ref_id_a = int32_t(id_map.try_insert(a.ref_name));
        contact_graph.try_insert_node(ref_id_a, 0);

        contact_graph.increment_coverage(ref_id_a, 1);

        for (size_t j=i+1; j<alignments.size(); j++) {
            auto& b = alignments[j];
            auto ref_id_b = int32_t(id_map.try_insert(b.ref_name));
            contact_graph.try_insert_node(ref_id_b, 0);
            contact_graph.try_insert_edge(ref_id_a, ref_id_b);
            contact_graph.increment_edge_weight(ref_id_a, ref_id_b, 1);
        }
    }
}


void parse_unpaired_bam_file(
        path bam_path,
        ContactGraph& contact_graph,
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
            update_contact_map(alignments, contact_graph, id_map);
            alignments.clear();
        }

        // No information about reference contig, this alignment is unusable
        if (a.ref_name.empty()){
            return;
        }

        // Optionally filter by the contig names. E.g. "PR" in shasta
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
            // Only allow reads with mapq > 0 and not secondary
            if (a.mapq >= min_mapq and a.is_primary()) {
                alignments.emplace_back(a);
            }
        }

        l++;
        prev_query_name = a.query_name;
    });
}


void phase_hic(path output_dir, path sam_path, path gfa_path, string required_prefix, int8_t min_mapq, size_t n_threads){
    Timer t;

    if (exists(output_dir)){
        throw runtime_error("ERROR: output directory exists already");
    }
    else {
        create_directories(output_dir);
    }

//    write_config(output_dir, sam_path, required_prefix, min_mapq, n_threads);

    // Id-to-name bimap for reference contigs
    IncrementalIdMap<string> id_map(false);

    GfaReader reader(gfa_path);

    // TODO: move this into domain of bdsg graph instead of GFA reader
    vector<Sequence> sequences;
    reader.for_each_sequence([&](string& name, string& sequence){
        id_map.try_insert(name);
        sequences.emplace_back(name, sequence);
    });

    // Datastructures to represent linkages from hiC
    ContactGraph contact_graph;
    vector <vector <int32_t> > adjacency;

    cerr << t << "Loading alignments as contact map..." << '\n';

    if (sam_path.extension() == ".bam"){
        parse_unpaired_bam_file(sam_path, contact_graph, id_map, required_prefix, min_mapq);
    }
    else{
        throw runtime_error("ERROR: unrecognized extension for BAM input file: " + sam_path.extension().string());
    }

//    HashGraph graph;

    for (auto& s: sequences){
        auto id = id_map.get_id(s.name);

        if (contact_graph.has_node(int32_t(id))) {
            contact_graph.set_node_length(int32_t(id), int64_t(s.size()));
        }
    }

    double sample_rate = 0.04;
    size_t k = 22;
    size_t n_iterations = 10;

    Hasher2 hasher(k, sample_rate, n_iterations, n_threads);

    hasher.hash(sequences);
    hasher.write_results(output_dir);

    // Before creating bubbles in the contact graph, write out all the node/edge data
    write_graph_data(output_dir, hasher, contact_graph, id_map);

    map<string,string> overlaps;

    hasher.get_symmetrical_matches(overlaps, 0.7);

    for (auto& [a,b]: overlaps){
        auto id_a = int32_t(id_map.get_id(a));
        auto id_b = int32_t(id_map.get_id(b));

        if (contact_graph.has_node(id_a) and contact_graph.has_node(id_b)){
            contact_graph.add_alt(id_a, id_b);
        }
    }

    path output_path = output_dir / "pairs.csv";
    ofstream file(output_path);

    if (not file.good() or not file.is_open()){
        throw runtime_error("ERROR: could not write file: " + output_path.string());
    }

    file << "Name" << ',' << "Match" << ',' << "Color" << '\n';
    for (auto& [a,b]: overlaps){
        file << a << ',' << b << ',' << "Cornflower Blue" << '\n';
        file << b << ',' << a << ',' << "Tomato" << '\n';
    }

    contact_graph.for_each_node([&](int32_t id, const Node& n){
        cerr << id << '\n';
        cerr << n << '\n';
    });

    vector <pair <int32_t,int8_t> > best_partitions;
    vector <int32_t> ids;
    atomic<double> best_score = std::numeric_limits<double>::min();
    atomic<size_t> job_index = 0;
    mutex phase_mutex;
    size_t m_iterations = 500;

    contact_graph.get_node_ids(ids);
    contact_graph.randomize_partitions();
    contact_graph.get_partitions(best_partitions);

    cerr << "Initial: " << std::flush;
    for (auto& [n,p]: best_partitions){
        cerr << '(' << id_map.get_name(n) << ',' << int(p) << ") ";
    }
    cerr << '\n';

    cerr << "start score: " << contact_graph.compute_total_consistency_score() << '\n';

    random_phase_search(contact_graph, ids, best_partitions, best_score, job_index, phase_mutex, m_iterations);

    contact_graph.set_partitions(best_partitions);

    cerr << "Final: " << std::flush;
    for (auto& [n,p]: best_partitions){
        cerr << '(' << id_map.get_name(n) << ',' << int(p) << ") ";
    }
    cerr << '\n';

    string suffix1 = "p" + required_prefix;
    string suffix2 = "m" + to_string(int(min_mapq));
    string suffix3 = "s" + to_string(best_score);

    cerr << t << "Writing phasing results to file... " << '\n';

    path contacts_output_path = output_dir / "contacts.csv";
    path phases_output_path = output_dir / "phases.csv";
    path config_output_path = output_dir / "config.csv";

    contact_graph.write_bandage_csv(phases_output_path, id_map);

    write_contact_map(contacts_output_path, contact_graph, id_map);

//    if (not gfa_path.empty()){
//        chain_phased_gfa(graph, id_map, bubble_graph, output_dir);
//    }

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

