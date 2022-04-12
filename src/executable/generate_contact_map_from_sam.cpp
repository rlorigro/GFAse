#include "Filesystem.hpp"
#include "CLI11.hpp"

using ghc::filesystem::path;
using CLI::App;

#include <unordered_map>
#include <unordered_set>
#include <stdexcept>
#include <fstream>
#include <utility>
#include <bitset>
#include <vector>
#include <array>
#include <set>

using std::unordered_set;
using std::unordered_map;
using std::runtime_error;
using std::to_string;
using std::ifstream;
using std::ofstream;
using std::string;
using std::vector;
using std::bitset;
using std::array;
using std::pair;
using std::stoi;
using std::cerr;
using std::set;


class SAMElement{
public:
    string read_name;
    string ref_name;
    int32_t line;
    int16_t flag;
    int8_t mapq;

    SAMElement();
    SAMElement(string& read_name, string& ref_name, int32_t line, int16_t flag, int8_t mapq);
    SAMElement(SAMElement&& other) noexcept;
    SAMElement(const SAMElement& other);
    bool is_first_mate() const;
    bool is_second_mate() const;
};


SAMElement::SAMElement(string& read_name, string& ref_name, int32_t line, int16_t flag, int8_t mapq):
    read_name(read_name),
    ref_name(ref_name),
    line(line),
    flag(flag),
    mapq(mapq)
{}


SAMElement::SAMElement(SAMElement&& other) noexcept:
    read_name(std::move(other.read_name)),
    ref_name(std::move(other.ref_name)),
    line(other.line),
    flag(other.flag),
    mapq(other.mapq)
{
    other.read_name.clear();
    other.ref_name.clear();
    other.line = 0;
    other.flag = 0;
    other.mapq = 0;
}


SAMElement::SAMElement(const SAMElement& other):
    read_name(other.read_name),
    ref_name(other.ref_name),
    line(other.line),
    flag(other.flag),
    mapq(other.mapq)
{}


SAMElement::SAMElement():
    read_name(),
    ref_name(),
    line(-1),
    flag(-1),
    mapq(-1)
{}


bool operator<(const SAMElement& a, const SAMElement& b){
    return a.line < b.line;
}


bool operator>(const SAMElement& a, const SAMElement& b){
    return a.line > b.line;
}


bool SAMElement::is_first_mate() const{
    return int16_t(flag) >> 6 & int16_t(1);
}


bool SAMElement::is_second_mate() const{
    return int16_t(flag) >> 7 & int16_t(1);
}

// Hash/compare using line number to guarantee unique
bool sam_comparator(const SAMElement& a, const SAMElement& b) {
    return a.line < b.line;
}


void print_mappings(const unordered_map <string, array <set <SAMElement, decltype(sam_comparator)*>, 2> >& mappings){
    for (const auto& [name,mates]: mappings){
        cerr << '\n';
        cerr << name << '\n';

        cerr << "First mates:" << '\n';
        for (auto& e: mates[0]){
            cerr << '\t' << e.ref_name << '\t' << bitset<sizeof(e.flag)*8>(e.flag) << '\n';
            cerr << '\t' << e.is_first_mate() << ' ' << e.is_second_mate() << '\n';
        }

        cerr << "Second mates:" << '\n';
        for (auto& e: mates[1]){
            cerr << '\t' << e.ref_name << '\t' << bitset<sizeof(e.flag)*8>(e.flag) << '\n';
            cerr << '\t' << e.is_first_mate() << ' ' << e.is_second_mate() << '\n';
        }
    }
}


void generate_contact_map_from_sam(path sam_path, string required_prefix, int8_t min_mapq){
    ifstream file(sam_path);

    if (not file.is_open() and file.good()){
        throw runtime_error("ERROR: could not read input file: " + sam_path.string());
    }

    char c;
    size_t n_delimiters = 0;
    int32_t n_lines = 0;

    string read_name;
    string flag_token;
    string ref_name;
    string mapq_token;

    unordered_map <string, array <set <SAMElement, decltype(sam_comparator)*>, 2> > mappings;

    while (file.get(c)){
        if (c == '\n'){
            n_delimiters = 0;

            bool valid_prefix = true;
            if (not required_prefix.empty()){
                for (size_t i=0; i<required_prefix.size(); i++){
                    if (ref_name[i] != required_prefix[i]){
                        valid_prefix = false;
                        break;
                    }
                }
            }

            if (valid_prefix) {
                auto mapq = int8_t(stoi(mapq_token));

                if (mapq >= min_mapq) {
                    SAMElement e(read_name, ref_name, n_lines, stoi(flag_token), mapq);
                    mappings[read_name][e.is_second_mate()].emplace(e);
                }
            }

//            cerr << n_lines << " " << read_name << " " << ref_name << '\n';

            read_name.clear();
            flag_token.clear();
            ref_name.clear();
            mapq_token.clear();

            n_lines++;
        }
        else if (c == '\t'){
            n_delimiters++;
        }
        else{
            if (n_delimiters == 0){
                read_name += c;
            }
            else if (n_delimiters == 1){
                flag_token += c;
            }
            else if (n_delimiters == 2){
                ref_name += c;
            }
            else if (n_delimiters == 4){
                mapq_token += c;
            }
        }
    }

    vector<string> to_be_deleted;
    for (auto& [name,mates]: mappings){

        bool has_first_mate = not mates[0].empty();
        bool has_second_mate = not mates[1].empty();

        if (not (has_first_mate and has_second_mate)){
            to_be_deleted.emplace_back(name);
        }
    }

    for (auto& item: to_be_deleted){
        mappings.erase(item);
    }

    unordered_map <string, unordered_map<string, int32_t> > contact_map;

    for (auto& [name,mates]: mappings) {
        size_t i = 0;

        for (auto& e: mates[0]){
            size_t j = 0;

            for (auto& e2: mates[1]){
                if (j >= i){
                    if (e.ref_name != e2.ref_name) {
                        contact_map[e.ref_name][e2.ref_name]++;
                        contact_map[e2.ref_name][e.ref_name]++;
                    }
                }
                j++;
            }
            i++;
        }
    }

    path output_path = sam_path;

    string suffix1 = "p" + required_prefix;
    string suffix2 = "m" + to_string(int(min_mapq));
    string suffix = suffix1 + "_" + suffix2 + "_contacts.csv";

    output_path.replace_extension(suffix);
    ofstream output_file(output_path);
    for (auto& [name,map2]: contact_map){
        for (auto& [name2,count]: map2){
            output_file << name << ',' << name2 << ',' << count << '\n';
        }
    }
}


int main (int argc, char* argv[]){
    path sam_path;
    string required_prefix;
    int8_t min_mapq = 0;

    CLI::App app{"App description"};

    app.add_option(
            "-i,--input",
            sam_path,
            "Path to SAM containing filtered, paired HiC reads")
            ->required();

    app.add_option(
            "-p,--prefix",
            required_prefix,
            "Prefix required in ref name for mapping to be counted");

    app.add_option(
            "-m,--min_mapq",
            min_mapq,
            "Minimum required mapq value for mapping to be counted");

    CLI11_PARSE(app, argc, argv);

    generate_contact_map_from_sam(sam_path, required_prefix, min_mapq);

    return 0;
}
