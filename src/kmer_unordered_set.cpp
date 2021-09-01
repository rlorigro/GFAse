#include <iostream>
#include <list>
#include <string>
#include "kmer_unordered_set.hpp"
#include <unordered_set>
#include <fstream>
// #include "Filesystem.hpp"

using namespace std;
// using ghc::filesystem::path;


void KmerSets::get_parent_kmer_sets(){
	// ask ryan about how to set up this file path for getting the file from the data folder
	// path script_path = __FILE__;
	// path project_directory = script_path.parent_path().parent_path().parent_path();

    // Get test region kmers
    string absolute_hg04_kmer_list_path = "/Users/mmmeredi/Desktop/paten/czi/mode2_unzipper/GFAse/data/hg04.all.homo.unique.kmer.1000.fa";
    // path absolute_kmerList = project_directory / relative_gfa_path;

	// Read from the text file
	ifstream MyReadFile(absolute_hg04_kmer_list_path);

	string line_txt;
	// Use a while loop to read the file line by line
	while (getline (MyReadFile, line_txt)) {
		// look for the > delimiting the kmer ID
		size_t found = line_txt.rfind('>');
		if (found!=string::npos){
			// insert the kmer into the respective kmer set
			// the line after the '>' in a fa
			getline (MyReadFile, line_txt);
			mom_kmer_set.insert(line_txt);
		}
	}

	// Close the file
	MyReadFile.close();
}

void KmerSets::fill_kmer_sets() {
	// hard code filled 
	mom_kmer_set.insert("GATTATTTTTACATTAAGGTTATCACCTCAAATCCTTTTTTAAAAATAGTCAGCA");
	mom_kmer_set.insert("GATTATTTTTACATTAAGGTTATCACCTCAAATCCTTTTTTAAAAATAGTCAGCC");
	mom_kmer_set.insert("GATTATTTTTACATTAAGGTTATCACCTCAAATCCTTTTTTAAAAATAGTCAGCT"); // child match
	
	dad_kmer_set.insert("AACTCTTCCCAGATGGATTTGAAACTTGAAATATGGAGGATAGAACTGTTACTGC");
	dad_kmer_set.insert("CACTCTTCCCAGATGGATTTGAAACTTGAAATATGGAGGATAGAACTGTTACTGC");
	dad_kmer_set.insert("GACTCTTCCCAGATGGATTTGAAACTTGAAATATGGAGGATAGAACTGTTACTGC"); // child match

}

// make a single kmer find
bool KmerSets::find_haplotype_single_kmer_count(string child_kmer) {
	mom_kmer_count = 0;
	dad_kmer_count = 0;

	if (mom_kmer_set.find(child_kmer) != mom_kmer_set.end()){
		mom_kmer_count++;
		}

	if (dad_kmer_set.find(child_kmer) != dad_kmer_set.end()){
		dad_kmer_count++;		}

	cout << " Checking single kmer: " << child_kmer << endl;
	cout << " found " << mom_kmer_count << " mom kmers." << endl ;
	cout << " found " << dad_kmer_count << " dad kmers." << endl << endl ;
	return 0;
}

bool KmerSets::find_haplotype_kmer_set_count(unordered_set <string> child_kmers) {
	unordered_set <string> :: iterator itr;
	mom_kmer_count = 0;
	dad_kmer_count = 0;
	for (itr = child_kmers.begin(); itr != child_kmers.end(); itr++)
		{
		if (mom_kmer_set.find(*itr) != mom_kmer_set.end()){
			mom_kmer_count++;
			}

		if (dad_kmer_set.find(*itr) != dad_kmer_set.end()){
			dad_kmer_count++;
			}
		}
	cout << "checking kmer set:" << endl;
	cout << " found " << mom_kmer_count << " mom kmers." << endl;
	cout << " found " << dad_kmer_count << " dad kmers." << endl << endl ;
	return 0;
}


