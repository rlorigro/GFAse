#include "bdsg/hash_graph.hpp"
#include "IncrementalIdMap.hpp"
#include "gfa_to_handle.hpp"
#include "GfaReader.hpp"
#include "KmerSets.hpp"
#include <unordered_set>
#include <string>

using namespace gfase;
using bdsg::HashGraph;

int main() {
	
	KmerSets<string> ks;
	ks.get_parent_kmer_sets();

	string path_string = "0-1";

	// test single kmer counting
	ks.increment_parental_kmer_count(path_string, "AAAAAAAAAAAAAAAAAAGGTGAAAGATCTGAACACCTCATTAATAAGATATACA"); // dad kmer
	ks.increment_parental_kmer_count(path_string, "AAAAAAAAAAAAAAAAAAGGTGCATGAAACATATGAAGCAAAAAGTGAAAGTCCC"); // dad kmer

	path_string = "0-0";
	ks.increment_parental_kmer_count(path_string, "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACAATTAAAAAA");
	
	path_string = "0-1";
	ks.increment_parental_kmer_count(path_string, "AAAAAAAAAAAAAAAAAAGGTGTCCATCCGAAAACCACCATTAAGAAACTCAGAC"); // dad kmer

	// test kmer set counting
	path_string = "1-1";
	unordered_set <string> child_kmer_set;
	child_kmer_set.insert("AAAAAAAAAAAAAAAAAAGGTGAAAGATCTGAACACCTCATTAATAAGATATACA"); //pat
	child_kmer_set.insert("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACACCCAAAAAA"); // mat
	child_kmer_set.insert("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACACAAAATAAA"); // mat
	
	ks.increment_parental_kmer_count(path_string, child_kmer_set);

	ks.print_component_parent_conf_matrix();

	// --- Use test data files to simulate running -- //
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();

    // Get test VCF path
    path relative_gfa_path = "data/simple_chain.gfa";
    path absolute_gfa_path = project_directory / relative_gfa_path;

    GfaReader reader(absolute_gfa_path);

    HashGraph graph;
    IncrementalIdMap<string> id_map;

    gfa_to_handle_graph(graph, id_map, absolute_gfa_path);


}
