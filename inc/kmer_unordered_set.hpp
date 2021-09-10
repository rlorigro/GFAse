#ifndef GFASE_KMERSETS_HPP
#define GFASE_KMERSETS_HPP

#include <iostream>
#include <list>
#include <string>
#include <map>
#include <unordered_set>
#include <fstream>
#include "Filesystem.hpp"

using namespace std;
using ghc::filesystem::path;

// Hashtable to implement a hashtable of kmers 

namespace gfase {


class KmerSets {
	/// Attributes ///
	private:
		// sets for each member of the trio
		unordered_set <string> hap1_kmer_set ; 
		unordered_set <string> hap2_kmer_set ; 
		// path to kmer parent files as input from command line
		path hap1_kmer_fa_path;
		path hap2_kmer_fa_path;
		// component and haplotype 
		int graph_component;
		int component_haplotype;
		static const int parent_hap1_int = 0; // hg03 paternal 
		static const int parent_hap2_int = 1; // hg04 maternal
		// < component,  [component_hap_path][parent_hap_int] >
		map<int, int [2][2]> component_map;

	/// Methods ///
	public:
		KmerSets();
		KmerSets(path hap1_kmer_fa_path_arg, path hap2_kmer_fa_path_arg);
		void load_file_into_unordered_set(path file_path, unordered_set <string>& set);
		void get_parent_kmer_sets();
		bool find_haplotype_kmer_set_count(string, unordered_set <string>);
		void parse_path_string(string);
		bool increment_parental_kmer_count(string, string );
		void print_component_parent_conf_matrix();
};

}

#endif //GFASE_KMERSETS_HPP
