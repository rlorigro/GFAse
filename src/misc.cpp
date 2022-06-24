#include "misc.hpp"
#include <map>

using std::map;

namespace gfase{


char get_reverse_complement(char c){
    if (c == 'A'){
        return 'T';
    }
    else if (c == 'C'){
        return 'G';
    }
    else if (c == 'G'){
        return 'C';
    }
    else if (c == 'T'){
        return 'A';
    }
    else {
        throw runtime_error("ERROR: uncomplementable character in sequence " + string(1,c));
    }
}


void get_reverse_complement(const string& fc, string& rc, size_t length){
    auto l = int64_t(length);

    for (int64_t i=l-1; i >= 0; i--){
        rc.append(1, get_reverse_complement(fc[i]));
    }
}



string join(const vector <string>& s, char delimiter){
    string joined_string;

    for (size_t i=0; i<s.size(); i++){
        if (i < s.size()-1){
            joined_string += s[i] + delimiter;
        }
        else{
            joined_string += s[i];
        }
    }

    return joined_string;
}


void run_command(const string& command){
    int exit_code = system(command.c_str());

    if (exit_code != 0){
        throw runtime_error("ERROR: command failed to run: " + command);
    }
}


path align(path output_dir, path ref_path, path query_path, size_t n_threads){
    path output_path;
    string prefix = query_path.filename().replace_extension("").string();
    string suffix = ref_path.filename().replace_extension("").string();

    replace(prefix.begin(), prefix.end(), '.', '_');
    replace(suffix.begin(), suffix.end(), '.', '_');

    output_path = prefix + "_VS_" + suffix + ".sam";
    output_path = absolute(output_dir) / output_path;

    // Set up arguments
    vector<string> arguments = {"minimap2",
                                "-a",
                                "-x", "asm20",
                                "-K", "10g",                    // New parameter for large batch size, better cpu %
                                "--eqx",
                                "-t", to_string(n_threads),
                                ref_path.string(),
                                query_path.string(),
                                "-o", output_path.string()
    };

    // Convert arguments to single string
    string command = join(arguments);
    cerr << "\nRUNNING: " << command << "\n";
    run_command(command);

    return output_path;
}


path sam_to_sorted_bam(path sam_path, size_t n_threads, bool remove_sam){

    vector<string> args = {
            "samtools", "sort",
            "-@", to_string(n_threads),
            "-o", "",
            ""
    };

    auto out_path = sam_path;
    out_path.replace_extension(".sorted.bam");
    args[5] = out_path;
    args[6] = sam_path;

    auto command = join(args);

    cerr << "\nRUNNING: " << command << "\n";
    run_command(command);

    if (remove_sam){
        remove(sam_path);
    }

    return out_path;
}


void get_query_lengths_from_fasta(path fasta_path, map<string,size_t>& query_lengths){
    ifstream file(fasta_path);

    auto iter = query_lengths.end();

    string name;
    string line;

    while (getline(file, line)){
        if (line[0] == '>'){
            // Trim any trailing tokens from the fasta header, keep only the name
            name = line.substr(1, line.find_first_of(" \t\n") - 1);

            auto r = query_lengths.emplace(name,0);

            if (r.second) {
                iter = r.first;
            }
            else{
                throw runtime_error("ERROR: failed to insert duplicate name into contig lengths: " + name);
            }
        }
        else {
            iter->second += line.size() - (isspace(line.back()));
        }
    }
}


void for_entry_in_csv(path csv_path, const function<void(const vector<string>& tokens, size_t line)>& f){
    ifstream file(csv_path);

    if (not file.is_open() or not file.good()){
        throw runtime_error("ERROR: could not read input file: " + csv_path.string());
    }

    char c;

    size_t n_delimiters = 0;
    size_t n_lines = 0;

    vector<string> empty_vector = {""};
    vector<string> tokens = empty_vector;

    while (file.get(c)){
        if (c == '\n'){
            n_delimiters = 0;

            if (tokens != empty_vector) {
                f(tokens, n_lines);
            }

            tokens = empty_vector;
            n_lines++;
        }
        else if (c == ','){
            n_delimiters++;
            tokens.emplace_back();
        }
        else {
            tokens.back() += c;
        }
    }

    if (c != '\n' and tokens != empty_vector){
        f(tokens, n_lines);
    }
}



}