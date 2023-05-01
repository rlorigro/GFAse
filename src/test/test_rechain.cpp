#include "MultiContactGraph.hpp"
#include "IncrementalIdMap.hpp"
#include "Sequence.hpp"
#include "gfa_to_handle.hpp"
#include "handle_to_gfa.hpp"
#include "hash_graph.hpp"
#include "HamiltonianChainer.hpp"
#include "Filesystem.hpp"

#include <iostream>
#include <cassert>

using namespace std;
using namespace gfase;
using namespace bdsg;

void load_alts_from_alignment_csv(path alignment_csv_path, MultiContactGraph& contact_graph, const IncrementalIdMap<string>& id_map){
    
    
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
                
                for (auto nid : {id_a, id_b}) {
                    if (!contact_graph.has_node(nid)) {
                        contact_graph.insert_node(nid);
                    }
                }
                
                bool symmetrical = stoi(symmetry_token);
                
                if (symmetrical and contact_graph.has_node(id_a) and contact_graph.has_node(id_b)){
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

void load_phasing_csv(path csv_path, const IncrementalIdMap<string>& id_map, MultiContactGraph& contact_graph)
{
    ifstream file(csv_path);
    
    if (not (file.is_open() and file.good())){
        throw runtime_error("ERROR: could not write to file: " + csv_path.string());
    }
    
    string line;
    vector<string> tokens = {""};
    
    size_t l = 0;
    while (getline(file, line)){
        if (l == 0){
            l++;
            continue;
        }
        
        for (auto& c: line){
            if (c == ','){
                tokens.emplace_back();
            }
            else{
                tokens.back() += c;
            }
        }
        
        auto& name = tokens[0];
        int phase = stoi(tokens[1]);
        auto node_id = id_map.get_id(name);
        contact_graph.set_partition(node_id, phase);
        
        
        tokens = {""};
        l++;
    }
}

int main(int argc, char** argv){
    
    assert(argc == 4);
    
    string graph_name = argv[1];
    string alignment_name = argv[2];
    string phasing_name = argv[3];
    
    path gfa_path = graph_name;
    path alignment_path = alignment_name;
    path phasing_path = phasing_name;
    
    IncrementalIdMap<string> id_map;
    HashGraph graph;
    Overlaps overlaps(graph);
    gfa_to_handle_graph(graph, id_map, overlaps, gfa_path, false, true);
    
    MultiContactGraph contact_graph;
    
    cerr << "loading alignments" << endl;
    load_alts_from_alignment_csv(alignment_path, contact_graph, id_map);

    cerr << "loading phasing" << endl;
    load_phasing_csv(phasing_path, id_map, contact_graph);
    
    HamiltonianChainer chainer;
    
    cerr << "chaining" << endl;
    
    chainer.generate_chain_paths(graph, id_map, contact_graph);
    cerr << "writing results" << endl;
    chainer.write_chaining_results_to_bandage_csv(path("./"), id_map, contact_graph);
    
    
    return 0;
}
