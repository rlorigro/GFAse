#include "Filesystem.hpp"
#include "CLI11.hpp"

using ghc::filesystem::path;
using CLI::App;

#include <unordered_map>
#include <vector>
#include <fstream>
#include <bitset>

using std::unordered_map;
using std::vector;
using std::ifstream;
using std::string;
using std::stoi;
using std::cerr;
using std::bitset;


class SAMElement{
public:
    string read_name;
    string ref_name;
    int16_t flag;
    int8_t mapq;

    SAMElement();
    SAMElement(string& read_name, string& ref_name, int16_t flag, int8_t mapq);
    bool is_first_mate();
    bool is_second_mate();
};


SAMElement::SAMElement(string& read_name, string& ref_name, int16_t flag, int8_t mapq):
    read_name(read_name),
    ref_name(ref_name),
    flag(flag),
    mapq(mapq)
{}


SAMElement::SAMElement():
    read_name(),
    ref_name(),
    flag(-1),
    mapq(-1)
{}


bool SAMElement::is_first_mate() {
    return int16_t(flag) >> 6 & int16_t(1);
}


bool SAMElement::is_second_mate() {
    return int16_t(flag) >> 7 & int16_t(1);
}


void generate_contact_map_from_sam(path sam_path){
    ifstream file(sam_path);

    char c;
    size_t n_delimiters = 0;

    string read_name;
    string flag_token;
    string ref_name;
    string mapq_token;

    unordered_map <string, vector<SAMElement> > mappings;

    while (file.get(c)){
        if (c == '\n'){
            n_delimiters = 0;

            SAMElement e(read_name, ref_name, stoi(flag_token), stoi(mapq_token));
            mappings[read_name].emplace_back(e);

            read_name.clear();
            flag_token.clear();
            ref_name.clear();
            mapq_token.clear();
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
    for (auto& item: mappings){

        if (item.second.size() == 1) {
            to_be_deleted.emplace_back(item.first);
        }
        else{
            bool has_first_mate = false;
            bool has_second_mate = false;

            for (auto& e: item.second){
                if (e.is_first_mate()){
                    has_first_mate = true;
                }
                if (e.is_second_mate()){
                    has_second_mate = true;
                }
            }

            if (not (has_first_mate and has_second_mate)){
                to_be_deleted.emplace_back(item.first);
            }
        }
    }

    for (auto& item: to_be_deleted){
        mappings.erase(item);
    }

    for (auto& item: mappings){
        cerr << '\n';
        cerr << item.first << '\n';
        for (auto& e: item.second){
            cerr << '\t' << e.ref_name << '\t' << bitset<sizeof(e.flag)*8>(e.flag) << '\n';
            cerr << e.is_first_mate() << ' ' << e.is_second_mate() << '\n';
        }
    }
}


int main (int argc, char* argv[]){
    path sam_path;

    CLI::App app{"App description"};

    app.add_option(
                    "-i,--input",
                    sam_path,
                    "Path to SAM containing filtered, paired HiC reads")
            ->required();

    CLI11_PARSE(app, argc, argv);

    generate_contact_map_from_sam(sam_path);

    return 0;
}
