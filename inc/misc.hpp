#ifndef GFASE_MISC_HPP
#define GFASE_MISC_HPP

#include "Filesystem.hpp"

using ghc::filesystem::create_directories;
using ghc::filesystem::exists;
using ghc::filesystem::path;

#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <utility>
#include <cstdio>
#include <limits>
#include <vector>
#include <array>
#include <set>
#include <map>


using std::numeric_limits;
using std::unordered_set;
using std::unordered_map;
using std::runtime_error;
using std::streamsize;
using std::exception;
using std::to_string;
using std::function;
using std::ifstream;
using std::ofstream;
using std::getline;
using std::remove;
using std::string;
using std::vector;
using std::array;
using std::pair;
using std::stoi;
using std::pair;
using std::hash;
using std::cerr;
using std::cref;
using std::ref;
using std::set;
using std::map;



namespace gfase{

string join(const vector <string>& s, char delimiter=' ');

void run_command(const string& argument_string);

path align(path output_dir, path ref_path, path query_path, size_t n_threads);

path sam_to_sorted_bam(path sam_path, size_t n_threads, bool remove_sam=true);

void get_query_lengths_from_fasta(path fasta_path, map<string,size_t>& query_lengths);

void for_entry_in_csv(path csv_path, const function<void(const vector<string>& tokens, size_t line)>& f);

char get_reverse_complement(char c);

void get_reverse_complement(const string& fc, string& rc, size_t length);


}


#endif //GFASE_MISC_HPP
