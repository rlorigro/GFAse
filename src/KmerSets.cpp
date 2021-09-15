#include "KmerSets.hpp"
#include "Filesystem.hpp"

#include <string>
#include <list>

using namespace std;
using ghc::filesystem::path;


namespace gfase {

pair<string, size_t> parse_path_string(string path_name, char delimiter){
    size_t index = path_name.find(delimiter);
    string component_name = path_name.substr(0,index);
    size_t component_haplotype = stoi(path_name.substr(index+1,path_name.length()));

    return {component_name, component_haplotype};
}

}
