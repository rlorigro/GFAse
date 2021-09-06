#ifndef GFASE_GFAREADER_HPP
#define GFASE_GFAREADER_HPP

#include <iostream>
#include <list>
#include <string>
#include <unordered_set>
#include <fstream>
#include "Filesystem.hpp"

using namespace std;
using ghc::filesystem::path;

// Hashtable to implement a hashtable of kmers 

class KmerSets {
	/// Attributes ///
	private:
		// sets for each member of the trio
		// unordered_set <string> childKmerSet ;
		unordered_set <string> mom_kmer_set ; // 482,003,862 hg04.all.homo.unique.ks.kmer.fa 	16G 
		unordered_set <string> dad_kmer_set ; // 716,947,716 hg03.all.homo.unique.kmer.fa 	    24G
		int mom_kmer_count;
		int dad_kmer_count;

	/// Methods ///
	public:
		void fill_kmer_sets();
		void load_file_into_unordered_set(path file_path, unordered_set <string>& set);
		void get_parent_kmer_sets();
		bool find_haplotype_kmer_set_count(unordered_set <string>);
		bool find_haplotype_single_kmer_count(string child_kmer);

};

#endif //GFASE_GFAREADER_HPP