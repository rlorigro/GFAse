#include "FixedBinarySequence.hpp"
#include "graph_utility.hpp"
#include "Filesystem.hpp"
#include "Phase.hpp"
#include "CLI11.hpp"

#include "bdsg/hash_graph.hpp"

#include <string>

using ghc::filesystem::path;

using gfase::FixedBinarySequence;
using gfase::HaplotypePathKmer;
using gfase::KmerSets;

using gfase::find_diploid_paths;
using gfase::parse_path_string;
using gfase::extend_paths;
using gfase::plot_graph;

using std::string;
using std::cout;
using std::cerr;


void construct_path_distance_map(const PathHandleGraph& graph, unordered_map<step_handle_t,size_t>& path_distance_map){
    graph.for_each_path_handle([&](const path_handle_t& p){
        size_t cumulative_length = 0;

        graph.for_each_step_in_path(p, [&](const step_handle_t& s){
            auto h = graph.get_handle_of_step(s);
            path_distance_map.emplace(s, cumulative_length);
            cumulative_length += graph.get_length(h);
        });
    });
}


void locate_kmer_matches(
        path gfa_path,
        size_t k,
        path paternal_kmers,
        path maternal_kmers,
        size_t min_path_length,
        set<string>& components,
        bool write_kmer_sequence,
        char path_delimiter = '.') {

    HashGraph graph;
    IncrementalIdMap<string> id_map;

    cerr << "Loading kmers into sets..." << '\n';
    KmerSets<FixedBinarySequence<uint64_t, 2> > ks(paternal_kmers, maternal_kmers);

    gfa_to_handle_graph(graph, id_map, gfa_path);

//    plot_graph(graph, "start_graph");

    cerr << "Identifying diploid paths..." << '\n';

    vector<path_handle_t> diploid_paths;

    if (not components.empty()){
        find_diploid_paths(graph, components, diploid_paths, '.');
    }
    else {
        find_diploid_paths(graph, diploid_paths);
    }

    cerr << "Extending paths by 1..." << '\n';

    extend_paths(graph);

    cerr << "Constructing path distance map..." << '\n';

    unordered_map<step_handle_t,size_t> path_distance_map;
    construct_path_distance_map(graph, path_distance_map);

    cerr << "Number of components in graph: " << graph.get_path_count() << '\n';

    path output_directory = gfa_path.parent_path() / (gfa_path.stem().string() + "_kmer_locations");

    cerr << "Writing kmer locations to: " << output_directory << '\n';
    create_directories(output_directory);

    path path_length_csv_path = output_directory / "path_lengths.csv";
    ofstream path_length_csv(path_length_csv_path);
    path_length_csv << "component_name" << ',' << "haplotype" << ',' << "length" << '\n';

    // Iterate paths and for each node, collect kmers if node is only covered by one path
    for (auto& p: diploid_paths) {
        string path_name = graph.get_path_name(p);

        HaplotypePathKmer kmer(graph, p, k);

        string component_name;
        size_t haplotype;

        tie(component_name, haplotype) = parse_path_string(path_name, path_delimiter);

        // Write out all the path lengths somewhere (for plotting convenience)
        uint64_t path_length = 0;
        graph.for_each_step_in_path(p, [&](const step_handle_t& s){
            path_length += graph.get_length(graph.get_handle_of_step(s));
        });

        if (path_length < min_path_length){
            // Skip to next path
            continue;
        }

        path_length_csv << component_name << ',' << haplotype << ',' << path_length << '\n';

        path output_subdirectory = output_directory / component_name;
        create_directories(output_subdirectory);

        path output_path = output_subdirectory / (to_string(haplotype) + ".csv");
        ofstream file(output_path);
        file << "path_index" << ',' << "is_paternal" << ',' << "is_maternal" << '\n';

        kmer.for_each_haploid_kmer([&](const deque<char>& sequence){
            FixedBinarySequence<uint64_t,2> s(sequence);

            // Get location of current kmer
            auto step = kmer.get_step_of_kmer_start();
            auto index = kmer.get_index_of_kmer_start();

            size_t path_position = path_distance_map.at(step) + index;

            // Get is_mat/is_pat
            bool is_paternal = ks.is_paternal(s);
            bool is_maternal = ks.is_maternal(s);

            file << path_position << ',' << int(is_paternal) << ',' << int(is_maternal);

            if (write_kmer_sequence){
                file << ',';
                for (const auto& c: sequence){
                    file << c;
                }
            }

            file << '\n';
        });
    }
}


int main (int argc, char* argv[]){
    path gfa_path;
    size_t k;
    size_t min_path_length;
    path paternal_kmers;
    path maternal_kmers;
    vector<string> c;
    set<string> components;
    bool write_kmer_sequence;

    CLI::App app{"App description"};

    app.add_option(
            "-i,--input_gfa",
            gfa_path,
            "Path to GFA containing phased non-overlapping segments")
            ->required();

    app.add_option(
            "-k,--kmer_size",
            k,
            "Length of kmer (k) to use")
            ->required();

    app.add_option(
            "-l,--min_length",
            min_path_length,
            "Minimum length of path to print information for")
            ->required();

    app.add_option(
            "-p,--paternal_kmers",
            paternal_kmers,
            "Paternal kmers in FASTA format")
            ->required();

    app.add_option(
            "-m,--maternal_kmers",
            maternal_kmers,
            "Maternal kmers in FASTA format")
            ->required();

    app.add_option(
            "-c,--components",
            c,
            "List of components to print (space separated)");

    app.add_flag(
            "-f,--full_kmer",
            write_kmer_sequence,
            "Boolean flag, add this option to write out the full kmer sequence for each kmer");

    CLI11_PARSE(app, argc, argv);

    for (const auto& item: c){
        components.emplace(item);
    }

    locate_kmer_matches(gfa_path, k, paternal_kmers, maternal_kmers, min_path_length, components, write_kmer_sequence);

    return 0;
}
