#ifndef GFASE_KMERSETS_HPP
#define GFASE_KMERSETS_HPP

#include <iostream>
#include <list>
#include <string>
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

		int hap1_kmer_count;
		int hap2_kmer_count;

	/// Methods ///
	public:
		KmerSets();
		KmerSets(path hap1_kmer_fa_path_arg, path hap2_kmer_fa_path_arg);
		void load_file_into_unordered_set(path file_path, unordered_set <string>& set);
		void get_parent_kmer_sets();
		bool find_haplotype_kmer_set_count(unordered_set <string>);
		pair<bool,bool> find_haplotype_single_kmer_count(string child_kmer);

};

}

#endif //GFASE_KMERSETS_HPP
