#ifndef GFASE_GFAREADER_HPP
#define GFASE_GFAREADER_HPP

#include "Filesystem.hpp"
#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <fstream>
#include <string>
#include <set>
#include <map>


using ghc::filesystem::path;
using std::unordered_map;
using std::unordered_set;
using std::function;
using std::ifstream;
using std::ofstream;
using std::string;
using std::vector;
using std::pair;
using std::set;
using std::map;


class GFAIndex{
public:
    /// Attributes ///
    char type;
    uint64_t offset;

    /// Methods ///
    GFAIndex(char type, uint64_t offset);
};


class GfaReader {
public:
    /// Attributes ///
    path gfa_path;
    path gfa_index_path;
    int gfa_file_descriptor;
    vector <GFAIndex> line_offsets;
    map <char, vector <size_t> > line_indexes_by_type;
    unordered_map <string, size_t> sequence_line_indexes_by_node;

    static const char EOF_CODE;

    /// Methods ///
    GfaReader(path gfa_path);
    ~GfaReader();
    void index();
    void read_index();
    void write_index_to_binary_file();
    void ensure_index_up_to_date();
    void map_sequences_by_node();
    void read_line(string& s, size_t index);
    uint64_t get_sequence_length(string node_name);
    void for_each_line_of_type(char type, const function<void(string& line)>& f);
    void for_each_sequence(const function<void(string& name, string& sequence)>& f);
    void for_each_link(const function<void(string& node_a, bool reversal_a, string& node_b, bool reversal_b, string& cigar)>& f);
    void for_each_path(const function<void(string& path_name, vector<string>& nodes, vector<bool>& reversals, vector<string>& cigars)>& f);
};


#endif //GFASE_GFAREADER_HPP
