#ifndef GFASE_KMERSETS_HPP
#define GFASE_KMERSETS_HPP

#include "graph_utility.hpp"
#include "Filesystem.hpp"
#include "spp.h"

#include <unordered_set>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <list>
#include <string>
#include <array>

using namespace std;
using spp::sparse_hash_set;
using ghc::filesystem::path;


namespace gfase {


char get_reverse_complement(char c);


void get_reverse_complement(const string& fc, string& rc, size_t length);


template <class T> class KmerSets {
	/// Attributes ///
	private:
		// Sets for each member of the trio
		sparse_hash_set <T> paternal_kmer_set;
		sparse_hash_set <T> maternal_kmer_set;

        // < component,  [component_hap_path][parent_hap_index] >
        unordered_map<string, array <array <double,2>, 2> > component_map;

        // Path to kmer parent files as input from command line
		path paternal_kmer_fa_path;
		path maternal_kmer_fa_path;
		double num_paternal_kmers;
		double num_maternal_kmers;

        size_t k = 0;

    public:
        static const size_t paternal_index = 0;
		static const size_t maternal_index = 1;

        // Character to use to split the path names into {component_name, haplotype}
        char path_delimiter;

    /// Methods ///
		KmerSets();

        // TODO: remove dependency on "path delimiter" for finding bubbles
		KmerSets(path paternal_kmer_fa_path_arg, path maternal_kmer_fa_path_args, char path_delimiter='.');
		float get_size_of_kmer_file(path file_path);
		void load_fasta_into_unordered_set(path file_path, sparse_hash_set<T>& s);
		void increment_parental_kmer_count(string path_hap_string, T child_kmer);
		void increment_parental_kmer_count(string path_name, unordered_set <T> child_kmers);
        void increment_parental_kmer_count(string component_name, size_t component_haplotype, T child_kmer);
        bool is_maternal(const T& kmer, const T& kmer_reverse_complement) const;
        bool is_paternal(const T& kmer, const T& kmer_reverse_complement) const;
        bool is_maternal(const T& kmer) const;
        bool is_paternal(const T& kmer) const;
		void normalize_kmer_counts();
		void print_component_parent_conf_matrix() const;
        void for_each_component_matrix(const function<void(const string& name, const array <array <double,2>, 2> component)>& f);
        double get_count(const string& component, size_t haplotype_index, size_t parental_index) const;
        size_t n_paternal_kmers() const;
        size_t n_maternal_kmers() const;
        void get_matrix(const string& component, array <array <double, 2>, 2>& matrix) const;
        size_t get_k();
};


template <class T> KmerSets<T>::KmerSets():
        num_paternal_kmers(0),
        num_maternal_kmers(0),
        path_delimiter('.')
{}


template <class T> KmerSets<T>::KmerSets(path paternal_kmer_fa_path, path maternal_kmer_fa_path, char path_delimiter):
        paternal_kmer_fa_path(paternal_kmer_fa_path),
        maternal_kmer_fa_path(maternal_kmer_fa_path),
        path_delimiter(path_delimiter)
{
    // Test files for hap1 and hap2 parents
    ifstream h1_test_stream(this->paternal_kmer_fa_path);
    if (not h1_test_stream.is_open()){
        throw runtime_error("ERROR: file could not be opened: " + this->paternal_kmer_fa_path.string());
    }

    ifstream h2_test_stream(this->maternal_kmer_fa_path);
    if (not h2_test_stream.is_open()){
        throw runtime_error("ERROR: file could not be opened: " + this->maternal_kmer_fa_path.string());
    }

    // Get the number of kmers for each parent
    num_paternal_kmers=get_size_of_kmer_file(paternal_kmer_fa_path);
    num_maternal_kmers=get_size_of_kmer_file(maternal_kmer_fa_path);
    cerr << " hap1 kmer file path: " << paternal_kmer_fa_path << "\n # kmers: " << num_paternal_kmers << endl;
    cerr << " hap2 kmer file path: " << maternal_kmer_fa_path << "\n # kmers: " << num_maternal_kmers << endl;

    // Fill the kmer sets
    load_fasta_into_unordered_set(paternal_kmer_fa_path, paternal_kmer_set);
    load_fasta_into_unordered_set(maternal_kmer_fa_path, maternal_kmer_set);
}


template <class T> size_t KmerSets<T>::get_k(){
    return k;
}


template <class T> size_t KmerSets<T>::n_paternal_kmers() const{
    return num_paternal_kmers;
}


template <class T> size_t KmerSets<T>::n_maternal_kmers() const{
    return num_maternal_kmers;
}


template <class T> float KmerSets<T>::get_size_of_kmer_file(path file_path){
    // # of kmers by (total file size/size of 2 lines)

    // Size of 2 fa lines = size of 1 kmer
    ifstream fileline(file_path, ios::binary );
    streampos fsize = 0;
    string line;

    // Move 1 kmer = 2 lines through file and store size
    getline(fileline, line);
    getline(fileline, line);
    fsize = fileline.tellg() - fsize;

    // Move to end
    fileline.seekg(0,ios_base::end);
    float size = fileline.tellg()/fsize;
    fileline.close();

    return size;
}


template <class T> void KmerSets<T>::load_fasta_into_unordered_set(path file_path, sparse_hash_set<T>& s){

    // Read from the text file
    ifstream KmerFile(file_path);

    string line;
    // Use a while loop to read the file line by line
    while (getline (KmerFile, line)) {
        if (line[0] == '>'){
            continue;
        }

        if (k == 0){
            k = line.size();
        }
        else if (line.size() != k){
            throw runtime_error("ERROR: kmer with unequal size found: " + line);
        }

        s.emplace(line);
    }
}


// Single kmer find
template <class T> void KmerSets<T>::increment_parental_kmer_count(string path_name, T child_kmer) {
    string graph_component;
    size_t component_haplotype;

    // This is done for each kmer - consider moving it out of this function
    tie(graph_component, component_haplotype) = parse_path_string(path_name, path_delimiter);

    increment_parental_kmer_count(graph_component, component_haplotype, child_kmer);
}


// Single kmer find
template <class T> void KmerSets<T>::increment_parental_kmer_count(
        string component_name,
        size_t component_haplotype,
        T child_kmer) {

    // Zero-initialize the arrays for each component
    auto iter = component_map.find(component_name);
    if (iter == component_map.end()){
        component_map.insert({component_name, {{{0, 0}, {0, 0}}}});
    }

    T child_kmer_rc;
    get_reverse_complement(child_kmer, child_kmer_rc, k);

    // Find child_kmer in each parent kmer_set
    component_map[component_name][component_haplotype][paternal_index] += is_paternal(child_kmer, child_kmer_rc);
    component_map[component_name][component_haplotype][maternal_index] += is_maternal(child_kmer, child_kmer_rc);
}


template <class T> bool KmerSets<T>::is_maternal(const T& kmer) const{
    T rc_kmer;
    get_reverse_complement(kmer, rc_kmer, k);

    bool found_forward = maternal_kmer_set.find(kmer) != maternal_kmer_set.end();
    bool found_reverse = maternal_kmer_set.find(rc_kmer) != maternal_kmer_set.end();

    return (found_forward or found_reverse);
}


template <class T> bool KmerSets<T>::is_paternal(const T& kmer) const{
    T rc_kmer;
    get_reverse_complement(kmer, rc_kmer, k);

    bool found_forward = paternal_kmer_set.find(kmer) != paternal_kmer_set.end();
    bool found_reverse = paternal_kmer_set.find(rc_kmer) != paternal_kmer_set.end();

    return (found_forward or found_reverse);
}


template <class T> bool KmerSets<T>::is_maternal(const T& kmer, const T& kmer_reverse_complement) const{
    bool found_forward = maternal_kmer_set.find(kmer) != maternal_kmer_set.end();
    bool found_reverse = maternal_kmer_set.find(kmer_reverse_complement) != maternal_kmer_set.end();

    return (found_forward or found_reverse);
}


template <class T> bool KmerSets<T>::is_paternal(const T& kmer, const T& kmer_reverse_complement) const{
    bool found_forward = paternal_kmer_set.find(kmer) != paternal_kmer_set.end();
    bool found_reverse = paternal_kmer_set.find(kmer_reverse_complement) != paternal_kmer_set.end();

    return (found_forward or found_reverse);
}


template <class T> void KmerSets<T>::increment_parental_kmer_count(string path_name, unordered_set <T> child_kmers) {
    for (auto& item: child_kmers){
        increment_parental_kmer_count(path_name,item);
    }
}


template <class T> void KmerSets<T>::normalize_kmer_counts(){
    // Divide by # of kmers for each parent

    cerr << "Normalize Component Map..\n";
    for (auto& [name, component]: component_map) {
        for (size_t j = 0; j < 2; j++) {
            component[j][paternal_index] /= num_paternal_kmers;
            component[j][maternal_index] /= num_maternal_kmers;
        }
    }
}


template <class T> double KmerSets<T>::get_count(const string& component, size_t haplotype_index, size_t parental_index) const{
    auto result = component_map.find(component);

    if (result == component_map.end()){
        cerr <<  "WARNING: component " + component + " not found in kmer counts, returning 0" << '\n';
        return 0;
    }

    return result->second[haplotype_index][paternal_index];
}


template <class T> void KmerSets<T>::get_matrix(const string& component, array <array <double, 2>, 2>& matrix) const{
    auto result = component_map.find(component);

    if (result == component_map.end()){
        cerr <<  "WARNING: component " + component + " not found in kmer counts, returning 0" << '\n';
        matrix[0][0] = 0;
        matrix[0][1] = 0;
        matrix[1][0] = 0;
        matrix[1][1] = 0;
    }
    else {
        matrix[0][0] = result->second[0][0];
        matrix[0][1] = result->second[0][1];
        matrix[1][0] = result->second[1][0];
        matrix[1][1] = result->second[1][1];
    }
}


template <class T> void KmerSets<T>::print_component_parent_conf_matrix() const{
    cout << "Number of graph components: " << component_map.size() << endl;

    for (const auto& [name, component]: component_map) {

        cout << "component: " << name << "\n  |pat\t|mat" << endl;

        for (size_t j = 0; j < 2; j++) {
            cout << j << " |" << component[j][paternal_index] << "\t|" << component[j][maternal_index] << "\t|" <<endl;
        }
        cout << endl;
    }
}

//, size_t hap, const size_t paternal_count, const size_t maternal_count, string end_delim
template <class T> void KmerSets<T>::for_each_component_matrix(const function<void(const string& name, const array <array <double,2>, 2> matrix)>& f){
    for (const auto& [name, matrix]: component_map) {
        // string end_delim = ",";
        // array of each component is hap0:[parent1,parent2], hap1:[parent1,parent2]
        // for (size_t j = 0; j < 2; j++) {
        //     if (j >0){
        //         end_delim = "\n";
        //     }
        f(name,matrix);
            //j,component[j][paternal_index],component[j][maternal_index],end_delim);
        }

}



    

}

#endif //GFASE_KMERSETS_HPP
