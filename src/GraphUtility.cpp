#include "GraphUtility.hpp"

namespace gfase {

void for_node_in_bfs(HandleGraph& graph, nid_t start_node, const function<void(const handle_t&)>& f) {
    unordered_set<nid_t> visited;
    queue<nid_t> q;

    q.emplace(start_node);
    visited.emplace(start_node);

    while (not q.empty()) {
        nid_t n = q.front();
        q.pop();

        auto h = graph.get_handle(n);
        f(h);

        graph.follow_edges(h, false, [&](const handle_t& other_handle) {
            auto other_node = graph.get_id(other_handle);

            // Attempt to add this to the set of visited nodes
            auto result = visited.emplace(other_node);

            // Check that this has NOT been visited before queuing it
            if (result.second) {
                q.emplace(other_node);
            }
        });

        graph.follow_edges(h, true, [&](const handle_t& other_handle) {
            auto other_node = graph.get_id(other_handle);

            // Attempt to add this to the set of visited nodes
            auto result = visited.emplace(other_node);

            // Check that this has NOT been visited before queuing it
            if (result.second) {
                q.emplace(other_node);
            }
        });
    }
}


void for_edge_in_bfs(HandleGraph& graph, nid_t start_node, const function<void(const handle_t& handle_a, const handle_t& handle_b)>& f) {
    unordered_set<nid_t> visited_nodes;
    unordered_map <handle_t, unordered_set<handle_t> > visited_edges;

    queue <handle_t> q;

    q.emplace(graph.get_handle(start_node));

    bool begin = true;

    while (not q.empty()) {
        handle_t h = q.front();
        q.pop();

        graph.follow_edges(h, false, [&](const handle_t& other_handle) {
            auto other_node = graph.get_id(other_handle);

            // Attempt to add this to the set of visited nodes
            auto node_result = visited_nodes.emplace(other_node);

            // Check that this has NOT been visited before queuing it
            if (node_result.second) {
                q.emplace(other_handle);
            }

            auto edge_result = visited_edges[h].insert(other_handle);

            if (edge_result.second){
                f(h, other_handle);
            }
        });

        graph.follow_edges(h, true, [&](const handle_t& other_handle) {
            auto other_node = graph.get_id(other_handle);

            // Attempt to add this to the set of visited nodes
            auto node_result = visited_nodes.emplace(other_node);

            // Check that this has NOT been visited before queuing it
            if (node_result.second) {
                q.emplace(other_handle);
            }

            auto edge_result = visited_edges[other_handle].emplace(h);

            if (edge_result.second){
                f(other_handle, h);
            }
        });

        if (begin){
            begin = false;
        }
    }
}


void for_each_connected_component(HandleGraph& graph, const function<void(unordered_set<nid_t>& connected_component)>& f) {
    unordered_set<nid_t> all_nodes;

    graph.for_each_handle([&](const handle_t& h) {
        auto n = graph.get_id(h);
        all_nodes.emplace(n);
    });

    while (not all_nodes.empty()) {
        unordered_set<nid_t> connected_component;
        auto iter = all_nodes.begin();

        for_node_in_bfs(graph, *iter, [&](const handle_t& h) {
            auto n = graph.get_id(h);

            connected_component.emplace(n);
            all_nodes.erase(n);
        });

        f(connected_component);
    }
}

/// For any 2 graphs with corresponding id maps, take a handle from the source graph and translate it into
/// an id in the destination graph, adding it to the destination id_map if it does not yet exist.
/// Then, return the id of the handle in source and destination as a pair
pair<nid_t,nid_t> translate_id(HandleGraph& source_graph,
                  IncrementalIdMap<string>& source_id_map,
                  HandleGraph& destination_graph,
                  IncrementalIdMap<string>& destination_id_map,
                  handle_t source_handle){

    auto id = source_graph.get_id(source_handle);
    auto name = source_id_map.get_name(id);

    cerr << "Translating: " << name << " with id " << id << '\n';

    nid_t translated_id;
    if (destination_id_map.exists(name)){
        translated_id = destination_id_map.get_id(name);
        cerr << "exists already: " << translated_id << '\n';
    }
    else{
        translated_id = destination_id_map.insert(name);
        cerr << "new id: " << translated_id << '\n';
    }

    return {id,translated_id};
}


void split_connected_components(
        PathHandleGraph& graph,
        IncrementalIdMap<string>& id_map,
        vector<HashGraph>& graphs,
        vector<IncrementalIdMap<string> >& id_maps) {

    unordered_set<nid_t> all_nodes;

    graph.for_each_handle([&](const handle_t& h) {
        auto n = graph.get_id(h);
        all_nodes.emplace(n);
    });

    while (not all_nodes.empty()) {
        cerr << "NEW CONNECTED COMPONENT" << '\n';
        graphs.emplace_back();
        id_maps.emplace_back();

        unordered_set<string> paths_to_be_copied;

        auto iter = all_nodes.begin();

        // Duplicate all the nodes
        for_node_in_bfs(graph, *iter, [&](const handle_t& h) {
            auto s = graph.get_sequence(h);
            nid_t id;
            nid_t other_id;
            tie(id, other_id) = translate_id(graph, id_map, graphs.back(), id_maps.back(), h);

            assert(not graphs.back().has_node(other_id));
            graphs.back().create_handle(s, other_id);

            if (id != *iter) {
                all_nodes.erase(id);
            }

            graph.for_each_step_on_handle(h, [&](const step_handle_t s){
                paths_to_be_copied.emplace(graph.get_path_name(graph.get_path_handle_of_step(s)));
            });
        });

        // Duplicate all the edges
        for_edge_in_bfs(graph, *iter, [&](const handle_t& handle_a, const handle_t& handle_b) {
            nid_t id_a;
            nid_t other_id_a;
            tie(id_a, other_id_a) = translate_id(graph, id_map, graphs.back(), id_maps.back(), handle_a);

            nid_t id_b;
            nid_t other_id_b;
            tie(id_b, other_id_b) = translate_id(graph, id_map, graphs.back(), id_maps.back(), handle_b);

            bool reversal_a = graph.get_is_reverse(handle_a);
            bool reversal_b = graph.get_is_reverse(handle_b);

            auto other_handle_a = graphs.back().get_handle(other_id_a, reversal_a);
            auto other_handle_b = graphs.back().get_handle(other_id_b, reversal_b);

            assert(not graphs.back().has_edge(other_handle_a, other_handle_b));
            graphs.back().create_edge(other_handle_a, other_handle_b);
        });

        // Duplicate all the paths
        for (auto& path_name: paths_to_be_copied){
            auto p = graph.get_path_handle(path_name);
            cerr << "copying path " << path_name << '\n';

            assert(not graphs.back().has_path(path_name));
            auto other_p = graphs.back().create_path_handle(path_name);

            graph.for_each_step_in_path(p, [&](const step_handle_t& s){
                auto h = graph.get_handle_of_step(s);

                nid_t id;
                nid_t other_id;
                tie(id, other_id) = translate_id(graph, id_map, graphs.back(), id_maps.back(), h);

                auto other_h = graph.get_handle(other_id, graph.get_is_reverse(h));

                graphs.back().append_step(other_p, other_h);
            });
        }

        all_nodes.erase(iter);
    }
}


void run_command(string& argument_string){
    int exit_code = system(argument_string.c_str());

    if (exit_code != 0){
        throw runtime_error("ERROR: command failed to run: " + argument_string);
    }
}


void plot_graph(const HandleGraph& graph, string filename_prefix){
    ofstream test_output(filename_prefix + ".gfa");
    handle_graph_to_gfa(graph, test_output);
    test_output.close();

    if (graph.get_node_count() < 200) {
        string command = "vg convert -g " + filename_prefix + ".gfa -p | vg view -d - | dot -Tpng -o "
                         + filename_prefix + ".png";

        cerr << "Running: " << command << '\n';

        run_command(command);
    }
}


void print_graph_paths(const PathHandleGraph& graph, const IncrementalIdMap<string>& id_map){
    graph.for_each_path_handle([&](const path_handle_t& p){
        auto path_name = graph.get_path_name(p);
        cerr << "Path " << path_name << '\n';

        graph.for_each_step_in_path(p, [&](const step_handle_t& s){
            auto h = graph.get_handle_of_step(s);
            auto id = graph.get_id(h);
            auto name = id_map.get_name(id);

            cerr << name << '\n';
        });
    });
}


}

