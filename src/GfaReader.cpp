#include "GfaReader.hpp"
#include "BinaryIO.hpp"
#include <iostream>
#include <fstream>
#include <string>

#include <sys/stat.h>
#include <ctime>
#include <cstdio>

using ::stat;
using std::difftime;
using std::stoi;
using std::cout;
using std::ofstream;
using std::runtime_error;


const char GfaReader::EOF_CODE = 'X';

GFAIndex::GFAIndex(char type, uint64_t offset){
    this->type = type;
    this->offset = offset;
}


void GfaReader::ensure_index_up_to_date(){
    struct stat index_file_stat;
    stat(gfa_index_path.c_str(), &index_file_stat);

    struct stat gfa_file_stat;
    stat(gfa_path.c_str(), &gfa_file_stat);

    cerr << "diff index_file_stat, gfa_file_stat: " << difftime(index_file_stat.st_ctime, gfa_file_stat.st_ctime) << '\n';

    if (difftime(index_file_stat.st_ctime, gfa_file_stat.st_ctime) < 0){
        throw runtime_error("ERROR index file is older than GFA: " + gfa_index_path.string());
    }
}


GfaReader::GfaReader(path gfa_path){
    this->gfa_path = gfa_path;
    this->gfa_index_path = gfa_path;
    this->gfa_index_path.replace_extension("gfai");
    this->gfa_file_descriptor = -1;

    // Test file
    ifstream test_stream(this->gfa_path);
    if (not test_stream.is_open()){
        throw runtime_error("ERROR: file could not be opened: " + this->gfa_path.string());
    }

    // Check if index exists, and generate one if necessary
    if (!exists(this->gfa_index_path)) {
        cerr << "No index found, generating .gfai for " << this->gfa_path << " ... ";

        this->index();
        cerr << "done\n";
    }
        // If index is found, load it
    else{
        cerr << "Found index, loading from disk: " << this->gfa_index_path << " ... ";

        this->read_index();
        cerr << "done\n";
    }

    this->map_sequences_by_node();
}


GfaReader::~GfaReader(){
    if (fcntl(this->gfa_file_descriptor, F_GETFD)){
        ::close(this->gfa_file_descriptor);
    }
}


void GfaReader::read_index(){
    // Open the input file.
    int file_descriptor = ::open(this->gfa_index_path.c_str(), O_RDONLY);

    // Verify it is working
    if(file_descriptor == -1) {
        throw runtime_error("ERROR: could not read " + this->gfa_index_path.string());
    }

    ensure_index_up_to_date();

    // Find file size in bytes, calculate number of entries in file
    off_t file_length = lseek(file_descriptor, 0, SEEK_END);
    off_t n_entries = file_length / (sizeof(uint64_t) + sizeof(char));

    // Initialize index to track position in file (AKA cursor)
    off_t byte_index = 0;

    // Read all the chunks and add them to the index data structures
    for (off_t i = 0; i < n_entries; i++){
        char c;
        uint64_t offset;

        pread_value_from_binary(file_descriptor,  c, byte_index);
        pread_value_from_binary(file_descriptor,  offset, byte_index);

        this->line_offsets.emplace_back(c, offset);
        this->line_indexes_by_type[c].emplace_back(this->line_offsets.size() - 1);
    }

    ::close(file_descriptor);
}


void GfaReader::write_index_to_binary_file(){
    ofstream index_file(this->gfa_index_path);

    for (auto& item: this->line_offsets){
        // Write a pair, denoting the line type and the byte offset for the start of that line
        write_value_to_binary(index_file, item.type);
        write_value_to_binary(index_file, item.offset);
    }
}


void GfaReader::index() {
    ifstream gfa_file(this->gfa_path);
    uint64_t byte_offset = 0;
    bool newline = true;
    char gfa_type_code;
    char c;

    // Find all the newlines in the GFA and store each line's byte offset in a vector. Additionally build a map which
    // lists all the positions in the index vector for each line type (e.g. S,L,H,U, etc.), so they can be iterated even
    // if they are not grouped or in order (which is not required by the GFA format spec)
    while (gfa_file.get(c)){
        if (c == '\n'){
            newline = true;
        }
        else{
            if (newline) {
                gfa_type_code = c;
                this->line_offsets.emplace_back(gfa_type_code, byte_offset);
                this->line_indexes_by_type[gfa_type_code].emplace_back(this->line_offsets.size() - 1);
                newline = false;
            }
        }
        byte_offset++;
    }

    // Append a placeholder to tell the total length of the file
    this->line_offsets.emplace_back(this->EOF_CODE, byte_offset);
    this->write_index_to_binary_file();
}


void GfaReader::read_line(string& s, size_t index){
    if (this->gfa_file_descriptor == -1){
        this->gfa_file_descriptor = ::open(this->gfa_path.c_str(), O_RDONLY);
    }

    off_t offset_start = this->line_offsets[index].offset;
    off_t offset_stop = this->line_offsets[index+1].offset;
    off_t length = offset_stop - offset_start;

    pread_string_from_binary(this->gfa_file_descriptor, s, length, offset_start);
}


void GfaReader::for_each_line_of_type(char type, const function<void(string& line)>& f) {
    if (this->line_indexes_by_type.count(type) == 0){
        return;
    }

    string line;
    for(auto& item: this->line_indexes_by_type.at(type)){
        read_line(line, item);
        f(line);
    }
}


void GfaReader::for_each_sequence(const function<void(string& name, string& sequence)>& f) {
    for_each_line_of_type('S', [&](string& line){
        size_t n_delimiters = 0;
        size_t position = 0;
        size_t start_position = 0;
        size_t stop_position = 0;
        string name;

        for (auto& c: line){
            if (isspace(c)){
                if (n_delimiters == 0){
                    start_position = position + 1;
                }
                else if (n_delimiters == 1){
                    stop_position = position;
                    name = line.substr(start_position, stop_position - start_position);

                    start_position = stop_position + 1;
                }
                else if (n_delimiters == 2){
                    stop_position = position;
                }

                n_delimiters++;
            }

            position++;
        }

        // Trim the line and send that in the "sequence" field to avoid copying the string again.
        line.erase(stop_position, line.size() - stop_position);
        line.erase(0, start_position);

        f(name, line);
    });
}


void GfaReader::for_each_link(const function<void(string& node_a, bool reversal_a, string& node_b, bool reversal_b, string& cigar)>& f) {
    size_t n_delimiters = 0;
    string node_a;
    string node_b;
    string cigar;
    bool reversal_a;
    bool reversal_b;

    for_each_line_of_type('L', [&](string& line){
        n_delimiters = 0;
        node_a.resize(0);
        node_b.resize(0);
        cigar.resize(0);

        for (auto& c: line){
            if (not isspace(c)){
                if (n_delimiters == 1){
                    node_a += c;
                }
                else if (n_delimiters == 2){
                    reversal_a = (c == '-');
                }
                else if (n_delimiters == 3){
                    node_b += c;
                }
                else if (n_delimiters == 4){
                    reversal_b = (c == '-');
                }
                else if (n_delimiters == 5){
                    cigar += c;
                }
            }
            else{
                n_delimiters++;
            }
        }

        f(node_a, reversal_a, node_b, reversal_b, cigar);
    });
}


void GfaReader::for_each_path(const function<void(string& path_name, vector<string>& nodes, vector<bool>& reversals, vector<string>& cigars)>& f) {
    size_t n_delimiters = 0;
    size_t n_subdelimiters = 0;
    string path_name;
    vector<string> nodes;
    vector<string> cigars;
    vector<bool> reversals;

    for_each_line_of_type('P', [&](string& line){
        n_delimiters = 0;
        n_subdelimiters = 0;
        path_name.resize(0);
        nodes.resize(0);
        nodes.emplace_back();
        cigars.resize(0);
        reversals.resize(0);
        size_t index = 0;

        for (auto& c: line){
            if (not isspace(c)){
//                cerr << c << ' ' << n_delimiters << ' ' << path_name << '\n';

                if (n_delimiters == 1){
                    path_name += c;
                }
                else if (n_delimiters == 2){
                    if (c != ','){
                        if ((c == '+' or c == '-')){
                            // TODO: check that reversal is always present
                            if (line[index+1] == ',') {
                                // Add the reversal for the current node
                                reversals.emplace_back((c == '-'));

                                // Add another entry to the back of the path to prepare for next node name
                                nodes.emplace_back();
                            }
                            if (isspace(line[index+1])) {
                                // Add the reversal for the current node ONLY
                                reversals.emplace_back((c == '-'));
                            }
                        }
                        else{
                            nodes.back() += c;
                        }
                    }
                    else{
                        n_subdelimiters++;
                    }
                }
                else if (n_delimiters == 3){
                    if (c == ','){
                        // Add another entry to the back of the path to prepare for next node name
                        cigars.emplace_back();
                    }
                    else{
                        if (cigars.empty()){
                            // "Cigars' is allowed to be empty if a path has only 1 node, so it should only be
                            // appended if a real cigar is found
                            cigars.emplace_back();
                        }

                        cigars.back() += c;
                    }
                }
            }
            else{
                n_delimiters++;
                n_subdelimiters = 0;
            }

            index++;
        }

        cerr << "----" << '\n';
        cerr << path_name << '\n';
        for (size_t i=0; i<nodes.size(); i++){
            cerr << nodes[i] << (reversals[i] ? '-' : '+') << ',';
        }
        cerr << '\n';

        for (auto& c: cigars){
            cerr << c << ',';
        }
        cerr << '\n';
        cerr << cigars.empty() << " " << cigars.size() << " " << nodes.size() << '\n';

        if (cigars.size() != nodes.size() - 1){
            throw runtime_error("ERROR: incorrect quantity of path cigars/overlaps for path: " + path_name);
        }

        f(path_name, nodes, reversals, cigars);
    });
}


void GfaReader::map_sequences_by_node(){
    cerr << "Mapping GFA S lines to node names... ";

    ifstream gfa_file(this->gfa_path);
    string token;
    char c = 0;

    // For every sequence line that has been indexed, jump to their offset in the file and read just the node names
    for (size_t i=0; i<this->line_indexes_by_type.at('S').size(); i++){
        auto line_index = this->line_indexes_by_type.at('S')[i];
        auto offset_start = this->line_offsets[line_index].offset + 2;  // Skip the line type character and tab
        token.resize(0);

        gfa_file.seekg(offset_start);

        // Parse the node name
        while (gfa_file.get(c)){
            if (c == '\t'){
                break;
            }
            token += c;
        }

        // Create the mapping that tells where in the file to find each node's sequence data
        this->sequence_line_indexes_by_node[token] = line_index;
        c = 0;
    }

    cerr << "done\n";
}


uint64_t GfaReader::get_sequence_length(string node_name){
    auto vector_index = this->sequence_line_indexes_by_node.at(node_name);
    auto start_index = this->line_offsets[vector_index].offset;

    char c;
    uint64_t length = 0;
    uint64_t n_separators = 0;
    ifstream file(this->gfa_path);

    file.seekg(start_index);

    // Count only characters in the sequence (not before or after)
    while(file.get(c)) {
        if (c == '\t' or c == '\n') {
            n_separators++;
        } else {
            if (n_separators == 2) {
                length++;
            }
            else if (n_separators > 2){
                break;
            }
        }
    }

    return length;
}

