#include "FixedBinarySequence.hpp"
#include "Phase.hpp"


namespace gfase{


void phase_haplotype_paths(path gfa_path, size_t k, path paternal_kmers, path maternal_kmers, char path_delimiter) {
    HashGraph graph;
    IncrementalIdMap<string> id_map;
    KmerSets <FixedBinarySequence <uint64_t, 2> > ks(paternal_kmers, maternal_kmers);

    gfa_to_handle_graph(graph, id_map, gfa_path);

    plot_graph(graph, "start_graph");

    cerr << "Identifying diploid paths..." << '\n';

    vector<string> diploid_path_names;
    find_diploid_paths(graph, diploid_path_names);

    cerr << "Extending paths by 1..." << '\n';

    vector<pair<path_handle_t, handle_t> > to_be_prepended;
    vector<pair<path_handle_t, handle_t> > to_be_appended;
    extend_paths(graph, to_be_prepended, to_be_appended);

    cerr << "Iterating path kmers..." << '\n';

    cerr << "\tNumber of components in graph: " << graph.get_path_count() << '\n';

    // Iterate paths and for each node, collect kmers if node is only covered by one path
    for (auto& path_name: diploid_path_names) {
        auto p = graph.get_path_handle(path_name);
        cerr << ">" << path_name << '\n';

        HaplotypePathKmer kmer(graph, p, k);

        string component_name;
        size_t haplotype;

        tie(component_name, haplotype) = parse_path_string(path_name, path_delimiter);

        kmer.for_each_haploid_kmer([&](const deque<char>& sequence){
            FixedBinarySequence<uint64_t,2> s(sequence);

            // Compare kmer to parental kmers
            ks.increment_parental_kmer_count(component_name, haplotype, s);
        });
    }

    cerr << "Un-extending paths..." << '\n';

    un_extend_paths(graph, to_be_prepended, to_be_appended);
    to_be_prepended.clear();
    to_be_appended.clear();

    // ks.normalize_kmer_counts();

    // Open file and print header
    ofstream component_matrix_outfile("component_matrix_outfile.csv");
    component_matrix_outfile << "component_name,hap_0_paternal_count,hap_0_maternal_count,hap_1_paternal_count,hap_1_maternal_count \n"; 
    
    ks.for_each_component_matrix([&](const string& name, const array <array <float,2>, 2> component){
        string end_delim = ",";
        for (size_t haplotype_index = 0; haplotype_index < 2; haplotype_index++) {
            if (haplotype_index == 0){
                component_matrix_outfile << name << ",";
            }
            else{
                end_delim = "\n";
            }
            component_matrix_outfile << component[haplotype_index][KmerSets<string>::paternal_index]
                                     << "," << component[haplotype_index][KmerSets<string>::maternal_index] << end_delim;
        }
    });

    cerr << "Unzipping..." << '\n';
    vector<HashGraph> connected_components;
    vector <IncrementalIdMap<string> > connected_component_ids;

    split_connected_components(graph, id_map, connected_components, connected_component_ids);

    for (size_t i=0; i<connected_components.size(); i++){
        cerr << "Component " << to_string(i) << '\n';
        print_graph_paths(connected_components[i], connected_component_ids[i]);

        unzip(connected_components[i], connected_component_ids[i]);

        string filename_prefix = "component_" + to_string(i) + "_unzipped";
        ofstream file(filename_prefix + ".gfa");
        handle_graph_to_gfa(connected_components[i], connected_component_ids[i], file);

        plot_graph(connected_components[i], filename_prefix);
    }

    ks.print_component_parent_conf_matrix();

    // TODO: for each path name in diploid paths, find the unzipped path name, and do bounded BFS on adjacent nodes
    // - Make no assumptions about what lies between diploid paths
    //    - Need to use DFS or DAG aligner with some scoring system for choosing path between components
    // - Filter out complex regions.. still need to decide on criteria for that
    //    - Depth of deepest snarl
    //    - Number of components touched by each unphased region must be < 3


}

}