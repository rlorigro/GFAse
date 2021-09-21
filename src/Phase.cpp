#include "FixedBinarySequence.hpp"
#include "Phase.hpp"


namespace gfase{


void phase_haplotype_paths(path gfa_path, size_t k, path paternal_kmers, path maternal_kmers, char path_delimiter) {
    HashGraph graph;
    IncrementalIdMap<string> id_map;
    KmerSets<FixedBinarySequence<uint64_t, 2> > ks(paternal_kmers, maternal_kmers);

    gfa_to_handle_graph(graph, id_map, gfa_path);

    // plot_graph(graph, "start_graph");

    cerr << "Identifying diploid paths..." << '\n';

    vector<path_handle_t> diploid_paths;
    find_diploid_paths(graph, diploid_paths);

    cerr << "Extending paths by 1..." << '\n';

    extend_paths(graph);

    cerr << "Iterating path kmers..." << '\n';

    cerr << "Number of components in graph: " << graph.get_path_count() << '\n';

    // Iterate paths and for each node, collect kmers if node is only covered by one path
    for (auto& p: diploid_paths) {
        string path_name = graph.get_path_name(p);
        cerr << ">" << path_name << '\n';

        HaplotypePathKmer kmer(graph, p, k);

        string component_name;
        size_t haplotype;

        tie(component_name, haplotype) = parse_path_string(path_name, path_delimiter);

        kmer.for_each_haploid_kmer([&](const deque<char>& sequence){
            FixedBinarySequence<uint64_t,2> s(sequence);

            // Compare kmer to parental kmers
            ks.increment_parental_kmer_count(component_name, haplotype, s);

//            for (const auto& c: sequence){
//                cerr << c;
//            }
//            cerr << '\n';
        });
    }

    // ks.normalize_kmer_counts();
    // ks.print_component_parent_conf_matrix();


    // open file and print header
    ofstream component_matrix_outfile;
    component_matrix_outfile.open("component_matrix_outfile.csv");
    component_matrix_outfile << "component_name,hap_0_paternal_count,hap_0_maternal_count,hap_1_paternal_count,hap_1_maternal_count \n"; 
    
    ks.for_each_component_matrix([&](const string& name, size_t hap, const size_t paternal_count, const size_t maternal_count, string end_delim){
            // print the name just once
            if (hap==0){
                component_matrix_outfile << name << ",";
            }
            component_matrix_outfile << paternal_count << "," << maternal_count << end_delim;
        });

    component_matrix_outfile.close();

}

}