#include "gfa_to_handle.hpp"

using bdsg::handle_t;
using bdsg::path_handle_t;
using bdsg::step_handle_t;
using bdsg::MutablePathHandleGraph;
using bdsg::HandleGraph;

namespace gfase {

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
        path gfa_file_path,
        bool ignore_singleton_paths){

    GfaReader gfa_reader(gfa_file_path);

    cerr << "Creating nodes..." << '\n';
    // Create all the nodes/sequences
    function<void(string&,string&)> f = [&](string& name, string& sequence) {
        // TODO: check if node name is empty or node sequence is empty
        auto id = parse_gfa_sequence_id(name, id_map);

        cerr << "graph: "<< &graph << endl;
        cerr << "id_map: "<< &id_map << endl;
        cerr << "name: "<< &name << endl;
        cerr << "id: "<< &id << endl;
        cerr << "sequence: "<< &sequence << endl;
        cerr << '\n';

        // Test graph is valid
        graph.get_node_count();
        graph.get_edge_count();
        graph.get_path_count();

        // Test id_map is valid
        for (auto& n: id_map.names){
            cerr << n << ',';
        }
        cerr << '\n';

        // Test sequence is valid
        for (auto& c: sequence){
            cerr << c;
        }
        cerr << '\n';

        // Test id is valid
        cerr << id << '\n';

        // Test name is valid
        cerr << name << '\n';

        graph.create_handle(sequence, id);
    };

    gfa_reader.for_each_sequence(f);

    cerr << "Creating edges..." << '\n';
    // Create all the edges between nodes
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
                                "contains non-empty overlap: " + cigar);
        }
    });

    cerr << "Creating paths..." << '\n';
    // Construct paths
    gfa_reader.for_each_path([&](string& path_name, vector<string>& nodes, vector<bool>& reversals, vector<string>& cigars){
        if (ignore_singleton_paths and nodes.size() == 1){
            return;
        }

        path_handle_t p = graph.create_path_handle(path_name);

        handle_t prev_handle;

        for (auto& cigar: cigars){
            if (not (cigar == "0M" or cigar == "*")){
                throw runtime_error("ERROR: cigar in path " + path_name + " contains non-empty overlap");
            }
        }

        for (size_t i=0; i<nodes.size(); i++){
            nid_t node_id;

            try {
                node_id = id_map.get_id(nodes[i]);
            }
            catch (const std::exception& e){
                cerr << e.what() << '\n';
                throw runtime_error("EEROR: node in path not found in GFA: " + nodes[i]);
            }

            handle_t handle = graph.get_handle(node_id, reversals[i]);

            if (i > 0 and not graph.has_edge(prev_handle, handle)){
                throw runtime_error("ERROR: graph has no edge between successive nodes in path: "
                                    + nodes[i-1] + (reversals[i-1] ? "-" : "+") + " -> " + nodes[i] + (reversals[i] ? "-" : "+"));
            }

            graph.append_step(p, handle);
            prev_handle = handle;
        }
    });
}


}