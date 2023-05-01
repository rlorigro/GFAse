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


void count_kmers(
        path gfa_path,
        size_t k,
        path paternal_kmers,
        path maternal_kmers,
        size_t min_path_length,
        char path_delimiter = '.') {

    HashGraph graph;
    IncrementalIdMap<string> id_map;
    Overlaps overlaps(graph);
    KmerSets <FixedBinarySequence <uint64_t, 2> > ks(paternal_kmers, maternal_kmers);

    gfa_to_handle_graph(graph, id_map, overlaps, gfa_path);

    cerr << "Identifying diploid paths..." << '\n';

    unordered_map<string, string> diploid_path_names;
    unordered_set<string> haploid_path_names;
    find_diploid_paths(graph, diploid_path_names, haploid_path_names);

    cerr << "Extending paths by 1..." << '\n';

    vector <pair<path_handle_t, handle_t> > to_be_prepended;
    vector <pair<path_handle_t, handle_t> > to_be_appended;
    extend_paths(graph, to_be_prepended, to_be_appended);

    cerr << "Iterating path kmers..." << '\n';

    cerr << "\tNumber of components in graph: " << graph.get_path_count() << '\n';

    // Iterate paths and for each node, collect kmers if node is only covered by one path
    for (auto& [path_name, other_path_name]: diploid_path_names) {
        auto p = graph.get_path_handle(path_name);

        uint64_t path_length = 0;
        graph.for_each_step_in_path(p, [&](const step_handle_t& s){
            path_length += graph.get_length(graph.get_handle_of_step(s));
        });

        if (path_length < min_path_length){
            // Skip to next path
            continue;
        }

        HaplotypePathKmer kmer(graph, p, k);

        string component_name;
        size_t haplotype;

        tie(component_name, haplotype) = parse_path_string(path_name, path_delimiter);

        cerr << path_name << " " << component_name << " " << haplotype << '\n';

        kmer.for_each_haploid_kmer([&](const deque<char>& sequence){
            FixedBinarySequence<uint64_t,2> s(sequence);

            // Compare kmer to parental kmers
            ks.increment_parental_kmer_count(component_name, haplotype, s);
        });
    }

    // Open file and print header
    ofstream component_matrix_outfile("kmer_counts.csv");
    component_matrix_outfile << "component_name,hap_0_paternal_count,hap_0_maternal_count,hap_1_paternal_count,hap_1_maternal_count \n";

    ks.for_each_component_matrix([&](const string& name, const array <array <double,2>, 2> matrix){
        component_matrix_outfile
                << name << ','
                << matrix[0][KmerSets<string>::paternal_index] << ','
                << matrix[0][KmerSets<string>::maternal_index] << ','
                << matrix[1][KmerSets<string>::paternal_index] << ','
                << matrix[1][KmerSets<string>::maternal_index] << '\n';
    });
}


int main (int argc, char* argv[]){
    path gfa_path;
    size_t k;
    size_t min_path_length;
    path paternal_kmers;
    path maternal_kmers;
    vector<string> c;

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

    CLI11_PARSE(app, argc, argv);

    count_kmers(gfa_path, k, paternal_kmers, maternal_kmers, min_path_length);

    return 0;
}
