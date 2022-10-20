#include "IncrementalIdMap.hpp"
#include "gfa_to_handle.hpp"
#include "graph_utility.hpp"
#include "MultiContactGraph.hpp"
#include "ContactGraph.hpp"
#include "Bipartition.hpp"
#include "hash_graph.hpp"
#include "Filesystem.hpp"
#include "Sequence.hpp"
#include "Hasher2.hpp"
#include "Timer.hpp"
#include "align.hpp"
#include "CLI11.hpp"
#include "Sam.hpp"
#include "Bam.hpp"
#include "minimap.h"

#include "SvgPlot.hpp"

using gfase::gfa_to_handle_graph;
using gfase::for_element_in_sam_file;


using gfase::random_multicontact_phase_search;
using gfase::construct_alignment_graph;
using gfase::unpaired_mappings_t;
using gfase::paired_mappings_t;
using gfase::contact_map_t;
using gfase::MultiContactGraph;
using gfase::ContactGraph;
using gfase::NonBipartiteEdgeException;

using gfase::AlignmentBlock;
using gfase::AlignmentChain;

using gfase::HashResult;
using gfase::Hasher2;

using gfase::IncrementalIdMap;
using gfase::Bipartition;
using gfase::SamElement;
using gfase::Sequence;
using gfase::Bubble;
using gfase::Timer;
using gfase::Node;
using gfase::Bam;

using bdsg::HashGraph;
using ghc::filesystem::path;
using CLI::App;


#include <unordered_map>
#include <thread>

using std::unordered_map;
using std::thread;
using std::cref;
using std::ref;


void load_alts_from_alignment_csv(path alignment_csv_path, MultiContactGraph& contact_graph, const IncrementalIdMap<string>& id_map){
    //     alignment_file << "name_a" << ',' << "name_b" << ',' << "total_matches" << ',' << "symmetrical" << ',' << "color" << '\n';


    ifstream file(alignment_csv_path);

    if (not (file.is_open() and file.good())){
        throw runtime_error("ERROR: could not read file: " + alignment_csv_path.string());
    }

    char c;

    string a;
    string b;
    string symmetry_token;

    size_t n_delimiters = 0;
    size_t n_lines = 0;

    while (file.get(c)){
        if (c == '\n'){
            if (n_lines > 0){
                auto id_a = int32_t(id_map.get_id(a));
                auto id_b = int32_t(id_map.get_id(b));

                bool symmetrical = stoi(symmetry_token);

                cerr << symmetry_token << ',' << symmetrical << ',' << a << ',' << b << '\n';

                if (symmetrical and contact_graph.has_node(id_a) and contact_graph.has_node(id_b)){
                    cerr << "adding alt: " << a << ',' << b << '\n';
                    contact_graph.add_alt(id_a, id_b);
                }
            }

            a.clear();
            b.clear();
            symmetry_token.clear();

            n_delimiters = 0;
            n_lines++;
        }
        else if (c == ','){
            n_delimiters++;
        }
        else{
            if (n_delimiters == 0){
                a += c;
            }
            else if (n_delimiters == 1){
                b += c;
            }
            else if (n_delimiters == 3){
                symmetry_token += c;
            }
            else if (n_delimiters > 4){
                throw runtime_error("ERROR: too many delimiters for line in file: " + alignment_csv_path.string());
            }
        }
    }
}


void rephase(
        path ids_csv_path,
        path alignment_csv_path,
        path contacts_path,
        path output_dir,
        size_t n_threads){

    if (exists(output_dir)){
        throw runtime_error("ERROR: output directory exists already");
    }
    else {
        create_directories(output_dir);
    }

    IncrementalIdMap<string> id_map(ids_csv_path);
    MultiContactGraph contact_graph(contacts_path, id_map);

    load_alts_from_alignment_csv(alignment_csv_path, contact_graph, id_map);

    phase_contacts(contact_graph, n_threads);

    path contacts_output_path = output_dir / "contacts.csv";
    path phases_output_path = output_dir / "phases.csv";

    contact_graph.write_bandage_csv(phases_output_path, id_map);
    contact_graph.write_contact_map(contacts_output_path, id_map);
}


int main (int argc, char* argv[]){
    path ids_csv_path;
    path alignment_csv_path;
    path contacts_path;
    path output_dir;
    size_t n_threads = 1;

    CLI::App app{"App description"};

    app.add_option(
            "-i,--ids",
            ids_csv_path,
            "Path to CSV describing node ids")
            ->required();

    app.add_option(
            "-a,--alignments",
            alignment_csv_path,
            "Path to CSV describing best node alignments")
            ->required();

    app.add_option(
            "-c,--contacts",
            contacts_path,
            "Path to CSV describing proximity linkage contacts")
            ->required();

    app.add_option(
            "-o,--output_dir",
            output_dir,
            "Path to (nonexistent) directory where output will be stored")
            ->required();

    app.add_option(
            "-t,--threads",
            n_threads,
            "Maximum number of threads to use");

    CLI11_PARSE(app, argc, argv);

    rephase(ids_csv_path, alignment_csv_path, contacts_path, output_dir, n_threads);

    return 0;
}


