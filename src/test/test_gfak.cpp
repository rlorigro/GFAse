#include <stdexcept>
#include <fstream>

#include "Filesystem.hpp"
#include "gfakluge.hpp"

using std::runtime_error;
using std::ifstream;
using std::cout;

using ghc::filesystem::path;
using ghc::filesystem::absolute;
using gfak::GFAKluge;


int main(){
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();

    // Get test GFA path
    path relative_gfa_path = "data/simple_chain.gfa";
    path absolute_gfa_path = project_directory / relative_gfa_path;

    cout << absolute_gfa_path << '\n';

    ifstream file(absolute_gfa_path);

    if (not (file.is_open() and file.good())){
        throw runtime_error("ERROR: file could not be read: " + absolute_gfa_path.string());
    }

    GFAKluge gfa_reader;

    gfa_reader.parse_gfa_file(file);

    for (auto& element: gfa_reader.get_name_to_seq()){
        cout << element.first << " " << element.second.sequence << '\n';
    }

    return 0;
}

