#include "handle_to_gfa.hpp"
#include "handlegraph/path_handle_graph.hpp"

using handlegraph::PathHandleGraph;
using handlegraph::path_handle_t;
using handlegraph::step_handle_t;

using std::runtime_error;

namespace gfase {


char get_reversal_character(const HandleGraph& graph, const handle_t& node){
    bool reversed = graph.get_is_reverse(node);

    if (reversed){
        return '-';
    }
    else{
        return '+';
    }
}


void write_node_to_gfa(const HandleGraph& graph, const handle_t& node, ostream& output_file){
    output_file << "S\t" << graph.get_id(node) << '\t' << graph.get_sequence(node) << '\n';
}


void write_node_to_gfa(const HandleGraph& graph, const IncrementalIdMap<string>& id_map, const handle_t& node, ostream& output_file){
    output_file << "S\t" << id_map.get_name(graph.get_id(node)) << '\t' << graph.get_sequence(node) << '\n';
}


void write_edge_to_gfa(const HandleGraph& graph, const edge_t& edge, ostream& output_file){
    output_file << "L\t" << graph.get_id(edge.first) << '\t' << get_reversal_character(graph, edge.first) << '\t'
                << graph.get_id(edge.second) << '\t' << get_reversal_character(graph, edge.second) << '\t'
                << "0M" << '\n';
}


void write_edge_to_gfa(const HandleGraph& graph, const IncrementalIdMap<string>& id_map, const edge_t& edge, ostream& output_file){
    output_file << "L\t" << id_map.get_name(graph.get_id(edge.first)) << '\t' << get_reversal_character(graph, edge.first) << '\t'
                << id_map.get_name(graph.get_id(edge.second)) << '\t' << get_reversal_character(graph, edge.second) << '\t'
                << "0M" << '\n';
}


void write_path_to_gfa(const PathHandleGraph& graph, const IncrementalIdMap<string>& id_map, const path_handle_t& path, ostream& output_file){
    string path_name = graph.get_path_name(path);
    size_t n_steps = graph.get_step_count(path);
    size_t i = 0;

    output_file << "P\t" << path_name << '\t';

    graph.for_each_step_in_path(path, [&](const step_handle_t& s){
        auto h = graph.get_handle_of_step(s);
        auto name = id_map.get_name(graph.get_id(h));

        output_file << name << (graph.get_is_reverse(h) ? '-' : '+');
        if (i < n_steps - 1){
            output_file << ',';
        }

        i++;
    });
    output_file << "\t";

    for (size_t j=0; j<n_steps-1; j++){
        output_file << "0M";
        if (j < n_steps - 2){
            output_file << ',';
        }
    }

    output_file << "\n";
}


/// With no consideration for directionality, just dump all the edges/nodes into GFA format
void handle_graph_to_gfa(const HandleGraph& graph, ostream& output_gfa){

    output_gfa << "H\tHVN:Z:1.0\n";

    graph.for_each_handle([&](const handle_t& node){
        write_node_to_gfa(graph, node, output_gfa);
    });

    graph.for_each_edge([&](const edge_t& edge){
        write_edge_to_gfa(graph, edge, output_gfa);
    });

    output_gfa.flush();
}


/// With no consideration for directionality, just dump all the edges/nodes into GFA format
void handle_graph_to_gfa(const PathHandleGraph& graph, const IncrementalIdMap<string>& id_map, ostream& output_gfa){

    output_gfa << "H\tHVN:Z:1.0\n";

    graph.for_each_handle([&](const handle_t& node){
        write_node_to_gfa(graph, id_map, node, output_gfa);
    });

    graph.for_each_edge([&](const edge_t& edge){
        write_edge_to_gfa(graph, id_map, edge, output_gfa);
    });

    graph.for_each_path_handle([&](const path_handle_t& path) {
        write_path_to_gfa(graph, id_map, path, output_gfa);
    });

    output_gfa << std::flush;
}

}
