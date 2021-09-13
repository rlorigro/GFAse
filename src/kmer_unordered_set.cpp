#include <iostream>
#include <list>
#include <string>
#include <array>
#include "kmer_unordered_set.hpp"
#include <unordered_set>
#include <fstream>
#include "Filesystem.hpp"

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

    // get the nimber of kmers for each parent 
    num_hap1_kmers=get_size_of_kmer_file(hap1_kmer_fa_path);
    num_hap2_kmers=get_size_of_kmer_file(hap2_kmer_fa_path);
    cerr << " hap1 kmer file path: " << hap1_kmer_fa_path << "\n # kmers: " << num_hap1_kmers << endl;
	cerr << " hap2 kmer file path: " << hap2_kmer_fa_path << "\n # kmers: " << num_hap2_kmers << endl;    
	// fill the kmer sets 
    load_file_into_unordered_set(hap1_kmer_fa_path, hap1_kmer_set);
    load_file_into_unordered_set(hap2_kmer_fa_path, hap2_kmer_set);

}



float KmerSets::get_size_of_kmer_file(path file_path){
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
	float size = fileline.tellg()/fsize;
	fileline.close();

	return size;
}

void KmerSets::load_file_into_unordered_set(path file_path, unordered_set <string>& set ){

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

    load_file_into_unordered_set(absolute_hap1_kmer_list_path,hap1_kmer_set);
    load_file_into_unordered_set(absolute_hap2_kmer_list_path,hap2_kmer_set);

}

void KmerSets::parse_path_string(string path_string){
	size_t string_delim = path_string.find('-');
	this->graph_component = stoi(path_string.substr(0,string_delim));
	this->component_haplotype = stoi(path_string.substr(string_delim+1,path_string.length()));
}

// single kmer find
bool KmerSets::increment_parental_kmer_count( string path_hap_string, string child_kmer) {

	// this is done for each kmer - consider moving it out of this function
	parse_path_string(path_hap_string);

	if (hap1_kmer_set.find(child_kmer) != hap1_kmer_set.end()){
		component_map[graph_component][component_haplotype][parent_hap1_index]++;
		}

	if (hap2_kmer_set.find(child_kmer) != hap2_kmer_set.end()){
		component_map[graph_component][component_haplotype][parent_hap2_index]++;	
		}

	return 0;
}

bool KmerSets::increment_parental_kmer_count(string path_name, unordered_set <string> child_kmers) {

	unordered_set <string> :: iterator itr;
	for (itr = child_kmers.begin(); itr != child_kmers.end(); itr++)
		{
			increment_parental_kmer_count(path_name,*itr);
		}

	return 0;
}

void KmerSets::normalize_kmer_counts(){
	// divide by # of kmers for each parent

	cerr << "Normalize Component Map..\n";
	for (size_t i = 0; i < component_map.size(); i++) {		
		for (size_t j = 0; j < 2; j++) {
			float raw_count = component_map[i][j][parent_hap1_index];
			component_map[i][j][parent_hap1_index] = (raw_count/num_hap1_kmers);

			raw_count = component_map[i][j][parent_hap2_index];
			component_map[i][j][parent_hap2_index] = (raw_count/num_hap2_kmers);

		}
	}
}

void KmerSets::print_component_parent_conf_matrix() {
	cout << "Number of graph components: " << component_map.size() << endl;
	
	for (size_t i = 0; i < component_map.size(); i++) {

		cout << "component: " << i << "\n  |pat\t|mat" << endl;
		
		for (size_t j = 0; j < 2; j++) {
			cout << j << " |" << component_map[i][j][parent_hap1_index] << "\t|" << component_map[i][j][parent_hap2_index] << "\t|" <<endl;
		}
		cout << endl;
	}
}


}
