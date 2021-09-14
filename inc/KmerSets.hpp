#ifndef GFASE_KMERSETS_HPP
#define GFASE_KMERSETS_HPP

#include "Filesystem.hpp"
#include "spp.h"

#include <iostream>
#include <list>
#include <string>
#include <array>
#include <map>
#include <unordered_set>
#include <fstream>

using namespace std;
using spp::sparse_hash_set;
using ghc::filesystem::path;


namespace gfase {


template <class T> class KmerSets {
	/// Attributes ///
	private:
		// Sets for each member of the trio
		sparse_hash_set <T> hap1_kmer_set ;
		sparse_hash_set <T> hap2_kmer_set ;

        // < component,  [component_hap_path][parent_hap_index] >
        map<size_t, array <array <float,2>, 2> > component_map;

        // Path to kmer parent files as input from command line
		path hap1_kmer_fa_path;
		path hap2_kmer_fa_path;
		float num_hap1_kmers;
		float num_hap2_kmers;

		static const size_t parent_hap1_index = 0;
		static const size_t parent_hap2_index = 1;

	/// Methods ///
	public:
		KmerSets();
		KmerSets(path hap1_kmer_fa_path_arg, path hap2_kmer_fa_path_args);
		float get_size_of_kmer_file(path file_path);
		void load_file_into_unordered_set(path file_path, unordered_set <T>& set);
		void get_parent_kmer_sets();
        pair<size_t, size_t> parse_path_string(string path_string);
		void increment_parental_kmer_count(string path_hap_string, T child_kmer);
		void increment_parental_kmer_count(string path_name, unordered_set <T> child_kmers);
		void normalize_kmer_counts();
		void print_component_parent_conf_matrix();
};


template <class T> KmerSets<T>::KmerSets():
        num_hap1_kmers(0),
        num_hap2_kmers(0)
{};


template <class T> KmerSets<T>::KmerSets(path hap1_kmer_fa_path, path hap2_kmer_fa_path):
        hap1_kmer_fa_path(hap1_kmer_fa_path),
        hap2_kmer_fa_path(hap2_kmer_fa_path)
{
    // Test files for hap1 and hap2 parents
    ifstream h1_test_stream(this->hap1_kmer_fa_path);
    if (not h1_test_stream.is_open()){
        throw runtime_error("ERROR: file could not be opened: " + this->hap1_kmer_fa_path.string());
    }

    ifstream h2_test_stream(this->hap2_kmer_fa_path);
    if (not h2_test_stream.is_open()){
        throw runtime_error("ERROR: file could not be opened: " + this->hap2_kmer_fa_path.string());
    }

    // Get the nimber of kmers for each parent
    num_hap1_kmers=get_size_of_kmer_file(hap1_kmer_fa_path);
    num_hap2_kmers=get_size_of_kmer_file(hap2_kmer_fa_path);
    cerr << " hap1 kmer file path: " << hap1_kmer_fa_path << "\n # kmers: " << num_hap1_kmers << endl;
    cerr << " hap2 kmer file path: " << hap2_kmer_fa_path << "\n # kmers: " << num_hap2_kmers << endl;

    // Fill the kmer sets
    load_file_into_unordered_set(hap1_kmer_fa_path, hap1_kmer_set);
    load_file_into_unordered_set(hap2_kmer_fa_path, hap2_kmer_set);
}


template <class T> float KmerSets<T>::get_size_of_kmer_file(path file_path){
    // # of kmers by (total file size/size of 2 lines)

    // Size of 2 fa lines = size of 1 kmer
    ifstream fileline(file_path, ios::binary );
    streampos fsize = 0;
    string fline;

    // Move 1 kmer = 2 lines through file and store size
    getline(fileline, fline);
    getline(fileline, fline);
    fsize = fileline.tellg() - fsize;

    // Move to end
    fileline.seekg(0,ios_base::end);
    float size = fileline.tellg()/fsize;
    fileline.close();

    return size;
}


template <class T> void KmerSets<T>::load_file_into_unordered_set(path file_path, unordered_set <T>& set ){

    // Read from the text file
    ifstream KmerFile(file_path);

    string line_txt;
    // Use a while loop to read the file line by line
    while (getline (KmerFile, line_txt)) {
        // Look for the > delimiting the kmer ID
        size_t found = line_txt.rfind('>');
        if (found!=string::npos){
            // Insert kmer (line after the '>') into set
            getline (KmerFile, line_txt);
            set.insert(line_txt);
        }
    }
    // Close the file
    KmerFile.close();
}


template <class T> void KmerSets<T>::get_parent_kmer_sets(){
    // Set up this file path for getting the file from the data folder
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


template <class T> pair<size_t, size_t> KmerSets<T>::parse_path_string(string path_name){
    size_t string_delim = path_name.find('-');
    size_t graph_component = stoi(path_name.substr(0,string_delim));
    size_t component_haplotype = stoi(path_name.substr(string_delim+1,path_name.length()));

    return {graph_component, component_haplotype};
}


// Single kmer find
template <class T> void KmerSets<T>::increment_parental_kmer_count( string path_name, T child_kmer) {
    size_t graph_component;
    size_t component_haplotype;

    // This is done for each kmer - consider moving it out of this function
    tie(graph_component, component_haplotype) = parse_path_string(path_name);

    // Zero-initialize the arrays for each component
    auto iter = component_map.find(graph_component);
    if (iter == component_map.end()){
        component_map.insert({graph_component, {{{0,0},{0,0}}}});
    }

    // Find child_kmer in each parent kmer_set
    if (hap1_kmer_set.find(child_kmer) != hap1_kmer_set.end()){
        component_map[graph_component][component_haplotype][parent_hap1_index]++;
    }

    if (hap2_kmer_set.find(child_kmer) != hap2_kmer_set.end()){
        component_map[graph_component][component_haplotype][parent_hap2_index]++;
    }
}


template <class T> void KmerSets<T>::increment_parental_kmer_count(string path_name, unordered_set <T> child_kmers) {
    for (auto& item: child_kmers){
        increment_parental_kmer_count(path_name,item);
    }
}


template <class T> void KmerSets<T>::normalize_kmer_counts(){
    // Divide by # of kmers for each parent

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


template <class T> void KmerSets<T>::print_component_parent_conf_matrix() {
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

#endif //GFASE_KMERSETS_HPP
