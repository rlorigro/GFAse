#include "kmer_unordered_set.hpp"
#include <unordered_set>
#include <string>

using namespace gfase;

int main() {
	
	KmerSets KS;
	KS.get_parent_kmer_sets();
	// change string into size_t 
	size_t pathi = 0;
	string path_string = "0-1"; 

	// test single kmer counting
	KS.increment_parental_kmer_count( path_string,"AAAAAAAAAAAAAAAAAAGGTGAAAGATCTGAACACCTCATTAATAAGATATACA"); // dad kmer
	KS.increment_parental_kmer_count( path_string,"AAAAAAAAAAAAAAAAAAGGTGCATGAAACATATGAAGCAAAAAGTGAAAGTCCC"); // dad kmer

	path_string = "0-0";
	KS.increment_parental_kmer_count( path_string,"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACAATTAAAAAA");
	
	path_string = "0-1";
	KS.increment_parental_kmer_count( path_string,"AAAAAAAAAAAAAAAAAAGGTGTCCATCCGAAAACCACCATTAAGAAACTCAGAC"); // dad kmer

	// test kmer set counting
	path_string = "1-1";
	unordered_set <string> child_kmer_set;
	child_kmer_set.insert("AAAAAAAAAAAAAAAAAAGGTGAAAGATCTGAACACCTCATTAATAAGATATACA"); //pat
	child_kmer_set.insert("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACACCCAAAAAA"); // mat
	child_kmer_set.insert("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACACAAAATAAA"); // mat
	
	KS.increment_parental_kmer_count(path_string, child_kmer_set);

	KS.print_component_parent_conf_matrix();

}
