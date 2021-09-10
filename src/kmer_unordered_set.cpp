#include <iostream>
#include <list>
#include <string>
#include "kmer_unordered_set.hpp"
#include <unordered_set>
#include <fstream>
#include "Filesystem.hpp"
// #include <bits/stdc++.h>

using namespace std;
using ghc::filesystem::path;

namespace gfase {

KmerSets::KmerSets(){};

KmerSets::KmerSets(path hap1_kmer_fa_path, path hap2_kmer_fa_path){
	// set kmer file paths
	this->hap1_kmer_fa_path = hap1_kmer_fa_path;
	this->hap2_kmer_fa_path = hap2_kmer_fa_path;

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
    load_file_into_unordered_set(hap1_kmer_fa_path, hap1_kmer_set, "hap1");
    load_file_into_unordered_set(hap2_kmer_fa_path, hap2_kmer_set, "hap2");
}

void KmerSets::get_size_of_kmer_file(path file_path, string hap){
	// # of kmers by (total file size/size of 2 lines)

	// size of 2 fa lines = size of 1 kmer
	ifstream fileline(file_path, ios::binary );
	streampos fsize = 0;
	string fline;
	// move 1 kmer = 2 lines through file and store size
	getline(fileline, fline);
	getline(fileline, fline);
    fsize = fileline.tellg() - fsize;
    // move to end
    fileline.seekg(0,ios_base::end);
	// set number of kmer variables
	if (hap=="hap1"){
		num_hap1_kmers=fileline.tellg()/fsize;
		cout << hap << " kmer file path: " << file_path << "\n # kmers: " << num_hap1_kmers << endl;
	}
	if (hap=="hap2"){
		num_hap2_kmers=fileline.tellg()/fsize;
		cout << hap << " kmer file path: " << file_path << "\n # kmers: " << num_hap2_kmers << endl;
	}

}

void KmerSets::load_file_into_unordered_set(path file_path, unordered_set <string>& set, string hap){
	
	get_size_of_kmer_file(file_path, hap);

	// Read from the text file
	ifstream KmerFile(file_path);

	string line_txt;
	// Use a while loop to read the file line by line
	while (getline (KmerFile, line_txt)) {
		// look for the > delimiting the kmer ID
		size_t found = line_txt.rfind('>');
		if (found!=string::npos){
			// insert kmer (line after the '>') into set
			getline (KmerFile, line_txt);
			set.insert(line_txt);
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

    load_file_into_unordered_set(absolute_hap1_kmer_list_path,hap1_kmer_set,"hap1");
    load_file_into_unordered_set(absolute_hap2_kmer_list_path,hap2_kmer_set,"hap2");

}


void KmerSets::print_component_parent_conf_matrix() {
	cout << "Number of graph components: " << component_map.size() << endl;
	
	for (int i = 0; i < component_map.size(); i++) {

		cout << "component: " << i << "\n   pat\tmat" << endl;
		
		for (int j = 0; j < 2; j++) {
			cout << j << " |" << component_map[i][j][parent_hap1_int] << "\t|" << component_map[i][j][parent_hap2_int] << "\t|" <<endl;
		}
		cout << endl;
	}
}

void KmerSets::parse_path_string(string path_string){
	int string_delim = path_string.find('-');
	this->graph_component = stoi(path_string.substr(0,string_delim));
	this->component_haplotype = stoi(path_string.substr(string_delim+1,path_string.length()));
}

// single kmer find
bool KmerSets::increment_parental_kmer_count( string path_hap_string, string child_kmer) {

	// this is done for each kmer - consider moving it out of this function
	parse_path_string(path_hap_string);

	if (hap1_kmer_set.find(child_kmer) != hap1_kmer_set.end()){
		component_map[graph_component][component_haplotype][parent_hap1_int]++;
		}

	if (hap2_kmer_set.find(child_kmer) != hap2_kmer_set.end()){
		component_map[graph_component][component_haplotype][parent_hap2_int]++;	
		}

	return 0;
}

void KmerSets::normalize_kmer_counts(){
	// divide by # of kmers for each parent

	cerr << "Normalize Component Map..\n";
	for (int i = 0; i < component_map.size(); i++) {		
		for (int j = 0; j < 2; j++) {
			float raw_count = component_map[i][j][parent_hap1_int];
			component_map[i][j][parent_hap1_int] = (raw_count/num_hap1_kmers);

			raw_count = component_map[i][j][parent_hap2_int];
			component_map[i][j][parent_hap2_int] = (raw_count/num_hap2_kmers);

		}
	}
}

bool KmerSets::find_haplotype_kmer_set_count(string path_hap_string, unordered_set <string> child_kmers) {

	unordered_set <string> :: iterator itr;
	for (itr = child_kmers.begin(); itr != child_kmers.end(); itr++)
		{
			increment_parental_kmer_count(path_hap_string,*itr);
		}

	return 0;
}


}
