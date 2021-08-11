#include "IncrementalIdMap.hpp"
#include "handle_to_gfa.hpp"
#include "Filesystem.hpp"
#include "GfaReader.hpp"
#include "CLI11.hpp"

#include "bdsg/hash_graph.hpp"

#include <string>

using gfase::IncrementalIdMap;
using gfase::handle_graph_to_gfa;

using ghc::filesystem::path;
using bdsg::HashGraph;
using bdsg::MutablePathMutableHandleGraph;
using bdsg::path_handle_t;
using bdsg::step_handle_t;
using bdsg::handle_t;

using std::string;
using std::cout;
using std::cerr;


nid_t parse_gfa_sequence_id(const string& s, IncrementalIdMap<string>& id_map) {
    nid_t id;

    // TODO rewrite this so it doesn't check for existing entry so many times...
    if (id_map.exists(s)){
        id = id_map.get_id(s);
    }
    else{
        id = id_map.insert(s);
    }

    return id;
}


void gfa_to_handle_graph(
        MutablePathMutableHandleGraph& graph,
        IncrementalIdMap<string>& id_map,
        unordered_map<nid_t,step_handle_t>& node_to_path_step,
        path gfa_file_path){

    GfaReader gfa_reader(gfa_file_path);

    gfa_reader.for_each_sequence([&](string& name, string& sequence){
        // TODO: check if node name is empty or node sequence is empty
        auto id = parse_gfa_sequence_id(name, id_map);
        graph.create_handle(sequence, id);
    });

    gfa_reader.for_each_link([&](string& node_a, bool reversal_a, string& node_b, bool reversal_b, string& cigar){
        const nid_t source_id = parse_gfa_sequence_id(node_a, id_map);
        const nid_t sink_id = parse_gfa_sequence_id(node_b, id_map);

        if (not graph.has_node(source_id)){
            throw runtime_error("ERROR: gfa link (" + node_a + "->" + node_b + ") "
                                "contains non-existent node: " + node_a);
        }

        if (not graph.has_node(sink_id)){
            throw runtime_error("ERROR: gfa link (" + node_a + "->" + node_b + ") "
                                "contains non-existent node: " + node_b);
        }

        // note: we're counting on implementations de-duplicating edges
        handle_t a = graph.get_handle(source_id, reversal_a);
        handle_t b = graph.get_handle(sink_id, reversal_b);
        graph.create_edge(a, b);

        if (cigar == "0M" or cigar == "*") {
            graph.create_edge(a, b);
        }
        else{
            throw runtime_error("ERROR: gfa link (" + node_a + "->" + node_b + ") "
                                "contains non empty overlap: " + cigar);
        }
    });

    gfa_reader.for_each_path([&](string& path_name, vector<string>& nodes, vector<bool>& reversals, vector<string>& cigars){
        path_handle_t p = graph.create_path_handle(path_name);

        for (size_t i=0; i<nodes.size(); i++){
            cout << nodes[i] << '\n';
            nid_t node_id = id_map.get_id(nodes[i]);
            handle_t handle = graph.get_handle(node_id, reversals[i]);
            step_handle_t step_handle = graph.append_step(p, handle);
            node_to_path_step.emplace(node_id, step_handle);

            //TODO: check if path is valid (edge exists)
        }

    });
}


void run_command(string& argument_string){
    int exit_code = system(argument_string.c_str());

    if (exit_code != 0){
        throw runtime_error("ERROR: command failed to run: " + argument_string);
    }
}


void unzip(path gfa_path){
    HashGraph graph;
    IncrementalIdMap<string> id_map;
    unordered_map<nid_t,step_handle_t> node_to_path_step;

    gfa_to_handle_graph(graph, id_map, node_to_path_step, gfa_path);

    // Output an image of the graph, can be uncommented for debugging
    {
        string test_path_prefix = "test_gfase_unedited";
        ofstream test_output(test_path_prefix + ".gfa");
        handle_graph_to_gfa(graph, test_output);
        test_output.close();

        if (graph.get_node_count() < 200) {
            string command = "vg convert -g " + test_path_prefix + ".gfa -p | vg view -d - | dot -Tpng -o "
                             + test_path_prefix + ".png";

            cerr << "Running: " << command << '\n';

            run_command(command);
        }
    }

    unordered_set<handle_t> to_be_destroyed;

    cout << graph.get_path_count() << '\n';
    graph.for_each_path_handle([&](const path_handle_t& p) {
        cout << graph.get_path_name(p) << "\n";

        string path_sequence;

        graph.for_each_step_in_path(p, [&](const step_handle_t s){
            handle_t h = graph.get_handle_of_step(s);
            nid_t n = graph.get_id(h);

            string sequence = graph.get_sequence(h);
            path_sequence += sequence;

            auto iter = node_to_path_step.find(n);
            if (iter != node_to_path_step.end()){
                node_to_path_step.erase(iter);
            }

            string name = id_map.get_name(n);
            cout << '\t' << name << " " << sequence << '\n';
        });

        // TODO: track provenance?
        handle_t haplotype_handle = graph.create_handle(path_sequence);
        auto path_start_handle = graph.get_handle_of_step(graph.path_begin(p));
        auto path_stop_handle = graph.get_handle_of_step(graph.path_back(p));

        graph.follow_edges(path_start_handle, true, [&](const handle_t& other){
            graph.create_edge(other, haplotype_handle);
        });

        graph.follow_edges(path_stop_handle, false, [&](const handle_t& other){
            graph.create_edge(haplotype_handle, other);
        });

        graph.for_each_step_in_path(p, [&](const step_handle_t s) {
            handle_t h = graph.get_handle_of_step(s);
            nid_t n = graph.get_id(h);

            auto iter = node_to_path_step.find(n);
            if (iter == node_to_path_step.end()) {
                to_be_destroyed.emplace(h);
            }
        });
    });

    // Output an image of the graph, can be uncommented for debugging
    {
        string test_path_prefix = "test_gfase_duplicated";
        ofstream test_output(test_path_prefix + ".gfa");
        handle_graph_to_gfa(graph, test_output);
        test_output.close();

        if (graph.get_node_count() < 200) {
            string command = "vg convert -g " + test_path_prefix + ".gfa -p | vg view -d - | dot -Tpng -o "
                             + test_path_prefix + ".png";

            cerr << "Running: " << command << '\n';

            run_command(command);
        }
    }

    for (auto& doomed_handle: to_be_destroyed){
        graph.destroy_handle(doomed_handle);
    }

    // Output an image of the graph, can be uncommented for debugging
    {
        string test_path_prefix = "test_gfase_unzipped";
        ofstream test_output(test_path_prefix + ".gfa");
        handle_graph_to_gfa(graph, test_output);
        test_output.close();

        if (graph.get_node_count() < 200) {
            string command = "vg convert -g " + test_path_prefix + ".gfa -p | vg view -d - | dot -Tpng -o "
                             + test_path_prefix + ".png";

            cerr << "Running: " << command << '\n';

            run_command(command);
        }
    }
}


int main (int argc, char* argv[]){
    path gfa_path;

    CLI::App app{"App description"};

    app.add_option("-i,--input_gfa", gfa_path, "Path to GFA containing phased non-overlapping segments")
            ->required();

    CLI11_PARSE(app, argc, argv);

    unzip(gfa_path);

    return 0;
}
