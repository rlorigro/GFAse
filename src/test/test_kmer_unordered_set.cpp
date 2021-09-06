#include "kmer_unordered_set.hpp"
#include <unordered_set>
#include <string>

int main() {
	unordered_set <string> child_kmer_set;
	child_kmer_set.insert("GATTATTTTTACATTAAGGTTATCACCTCAAATCCTTTTTTAAAAATAGTCAGCA");
	child_kmer_set.insert("GGACTCTTCCCAGATGGATTTGAAACTTGAAATATGGAGGATAGAACTGTTACTG");
	child_kmer_set.insert("GACTCTTCCCAGATGGATTTGAAACTTGAAATATGGAGGATAGAACTGTTACTGC");

	KmerSets KS;
	KS.get_parent_kmer_sets();
	// KS.fill_kmer_sets();  // with kmers that match the child set here
	KS.find_haplotype_single_kmer_count("GACTCTTCCCAGATGGATTTGAAACTTGAAATATGGAGGATAGAACTGTTACTGC"); // dad kmer
	KS.find_haplotype_single_kmer_count("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATCG"); // mom kmer
	KS.find_haplotype_kmer_set_count(child_kmer_set);

}
