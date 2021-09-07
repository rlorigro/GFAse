#include <iostream>
#include <list>
#include <string>
#include "kmer_unordered_set.hpp"
#include <unordered_set>
#include <fstream>
#include "Filesystem.hpp"

using namespace std;
using ghc::filesystem::path;

KmerSets::KmerSets(){};

KmerSets::KmerSets(path hap1_kmer_fa_path_arg, path hap2_kmer_fa_path_arg){
	// set kmer file paths
	hap1_kmer_fa_path = hap1_kmer_fa_path_arg;
	hap2_kmer_fa_path = hap2_kmer_fa_path_arg;

	// Test files for hap1 and hap2 parents
    ifstream h1_test_stream(this->hap1_kmer_fa_path);
    if (not h1_test_stream.is_open()){
        throw runtime_error("ERROR: file could not be opened: " + this->hap1_kmer_fa_path.string());
    }

    ifstream h2_test_stream(this->hap2_kmer_fa_path);
    if (not h2_test_stream.is_open()){
        throw runtime_error("ERROR: file could not be opened: " + this->hap2_kmer_fa_path.string());
    }

    // fill the kmer sets 
    load_file_into_unordered_set(hap1_kmer_fa_path, hap1_kmer_set);
    load_file_into_unordered_set(hap2_kmer_fa_path, hap2_kmer_set);
}

void KmerSets::load_file_into_unordered_set(path file_path, unordered_set <string>& set){
	// Read from the text file
	ifstream KmerFile(file_path);
	cout << "kmer path: " << file_path << endl << endl;

	string line_txt;
	// Use a while loop to read the file line by line
	while (getline (KmerFile, line_txt)) {
		// look for the > delimiting the kmer ID
		size_t found = line_txt.rfind('>');
		if (found!=string::npos){
			// insert the kmer into the respective kmer set
			// the line after the '>' in a fa
			getline (KmerFile, line_txt);
			set.insert(line_txt);
			// cout << " kmer " << line_txt << endl << endl;
		}
	}
	// Close the file
	KmerFile.close();
	
}

void KmerSets::get_parent_kmer_sets(){
	// set up this file path for getting the file from the data folder
	path script_path = __FILE__;
	path project_directory = script_path.parent_path().parent_path(); // this path is different than the one ryan uses in GfaReader and I'm not sure why

    // Get test parent1 kmers
    path relative_hap1_kmer_list_path = "data/hg03.all.homo.unique.kmer.1000.fa";
    path absolute_hap1_kmer_list_path = project_directory / relative_hap1_kmer_list_path;

    // Get test parent2 kmers
    path relative_hap2_kmer_list_path = "data/hg04.all.homo.unique.kmer.1000.fa";
    path absolute_hap2_kmer_list_path = project_directory / relative_hap2_kmer_list_path;

    load_file_into_unordered_set(absolute_hap1_kmer_list_path,hap1_kmer_set);
    load_file_into_unordered_set(absolute_hap2_kmer_list_path,hap2_kmer_set);

}

// single kmer find
bool KmerSets::find_haplotype_single_kmer_count(string child_kmer) {
	hap1_kmer_count = 0;
	hap2_kmer_count = 0;

	if (hap1_kmer_set.find(child_kmer) != hap1_kmer_set.end()){
		hap1_kmer_count++;
		}

	if (hap2_kmer_set.find(child_kmer) != hap2_kmer_set.end()){
		hap2_kmer_count++;		
		}

	cout << " Checking single kmer: " << child_kmer << endl;
	cout << " found " << hap1_kmer_count << " hap1 (hg03) kmers." << endl ;
	cout << " found " << hap2_kmer_count << " hap2 (hg04) kmers." << endl << endl ;
	return 0;
}

// kmer set find
bool KmerSets::find_haplotype_kmer_set_count(unordered_set <string> child_kmers) {
	unordered_set <string> :: iterator itr;
	hap1_kmer_count = 0;
	hap2_kmer_count = 0;

	for (itr = child_kmers.begin(); itr != child_kmers.end(); itr++)
		{
		if (hap1_kmer_set.find(*itr) != hap1_kmer_set.end()){
			hap1_kmer_count++;
			}

		if (hap2_kmer_set.find(*itr) != hap2_kmer_set.end()){
			hap2_kmer_count++;
			}
		}
	cout << "checking kmer set." << endl;
	cout << " found " << hap1_kmer_count << " hap1 (hg03) kmers." << endl;
	cout << " found " << hap2_kmer_count << " hap2 (hg04) kmers." << endl << endl ;
	return 0;
}


