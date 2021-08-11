#include <GfaReader.hpp>
#include <iostream>

using std::cerr;
using std::cerr;
using std::stringstream;

void test_gfa(GfaReader& reader){
    for (auto& [gfa_type_code, indexes]: reader.line_indexes_by_type){
        for (auto& i: indexes) {
            cerr << gfa_type_code << '\t' << i << '\t' << reader.line_offsets[i].offset << '\n';
        }
    }

    cerr << reader.line_offsets.back().type << '\t' << reader.line_offsets.back().offset << '\n';

    cerr << reader.sequence_line_indexes_by_node.size() << '\n';
    for (auto& item: reader.sequence_line_indexes_by_node){
        cerr << item.first << " " << item.second << '\n';
    }

    cerr << "TESTING sequence iterator\n";
    reader.for_each_sequence([&](string& name, string& sequence){
        cerr << name << ' ' << sequence << '\n';
    });

    cerr << "TESTING link iterator\n";
    reader.for_each_link([&](string& node_a, bool reversal_a, string& node_b, bool reversal_b, string& cigar){
        cerr << node_a << " " << reversal_a << " " << node_b << " " << reversal_b << " " << cigar << '\n';
    });

    cerr << "TESTING path iterator\n";
    reader.for_each_path([&](string& path_name, vector<string>& nodes, vector<bool>& reversals, vector<string>& cigars){
        cerr << path_name << '\n';

        for (auto& node: nodes){
            cerr << node << ',';
        }
        cerr << '\n';
        for (auto reversal: reversals){
            cerr << (reversal ? '-' : '+') << ',';
        }
        cerr << '\n';
        for (auto& cigar: cigars){
            cerr << cigar << ',';
        }
        cerr << '\n';
    });

}


int main(){
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();

    // Get test VCF path
    path relative_gfa_path = "data/test_gfa1.gfa";
    path absolute_gfa_path = project_directory / relative_gfa_path;

    GfaReader reader(absolute_gfa_path);
    test_gfa(reader);

    // Get test VCF path
    path relative_gfa_path_2 = "data/simple_chain.gfa";
    path absolute_gfa_path_2 = project_directory / relative_gfa_path_2;

    GfaReader reader_2(absolute_gfa_path_2);
    test_gfa(reader_2);

    return 0;
}

