#ifndef GFASE_KMERSETS_HPP
#define GFASE_KMERSETS_HPP

#include <iostream>
#include <list>
#include <string>
#include <array>
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
		float num_hap1_kmers;
		float num_hap2_kmers;
		// component and haplotype 
		size_t graph_component;
		size_t component_haplotype;
		static const size_t parent_hap1_index = 0; // hg03 paternal 
		static const size_t parent_hap2_index = 1; // hg04 maternal
		// < component,  [component_hap_path][parent_hap_index] >
		map<size_t, array <array <float,2>, 2>> component_map;

	/// Methods ///
	public:
		KmerSets();
		KmerSets(path hap1_kmer_fa_path_arg, path hap2_kmer_fa_path_args);
		void load_file_into_unordered_set(path , unordered_set <string>& );
		float get_size_of_kmer_file(path);
		void get_parent_kmer_sets();
		bool parental_kmer_count_for_kmer_set(string, unordered_set <string>);
		void parse_path_string(string);
		bool increment_parental_kmer_count(string, string );
		void normalize_kmer_counts();
		void print_component_parent_conf_matrix();
};

}

#endif //GFASE_KMERSETS_HPP
