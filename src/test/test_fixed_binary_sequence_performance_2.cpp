#include "FixedBinarySequence.hpp"
#include "GfaReader.hpp"
#include "Filesystem.hpp"
#include "CLI11.hpp"

using ghc::filesystem::path;
using gfase::FixedBinarySequence;

#include <iostream>

using std::cerr;
using std::cerr;
using std::stringstream;


int test(path gfa_path){
    GfaReader reader(gfa_path);

    size_t k = 12;
    FixedBinarySequence<int32_t,1> s;


    reader.for_each_sequence([&](string& name, string& sequence){
        cerr << name << '\n';
        size_t total = 0;
        string out;
        for (auto& c: sequence){
            s.shift(c, k);

            total += s.sequence.front();
        }
        cerr << sequence.size() << ' ' << total << '\n';
    });

    return 0;
}


int main (int argc, char* argv[]){
    path gfa_path;

    CLI::App app{"App description"};

    app.add_option(
            "-i,--input",
            gfa_path,
            "Path to GFA")
            ->required();

    CLI11_PARSE(app, argc, argv);

    test(gfa_path);

    return 0;
}
