#include "graph_utility.hpp"

namespace gfase {


pair<string, size_t> parse_path_string(string path_name, char delimiter){
    size_t index = path_name.find(delimiter);
    string component_name = path_name.substr(0,index);
    size_t component_haplotype = stoi(path_name.substr(index+1,path_name.length()));

    return {component_name, component_haplotype};
}


void for_node_in_bfs(const HandleGraph& graph, nid_t start_node, const function<void(const handle_t&)>& f) {
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


void for_edge_in_bfs(const HandleGraph& graph, nid_t start_node, const function<void(const handle_t& handle_a, const handle_t& handle_b)>& f) {
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
pair<nid_t,nid_t> translate_id(const HandleGraph& source_graph,
                  const IncrementalIdMap<string>& source_id_map,
                  HandleGraph& destination_graph,
                  IncrementalIdMap<string>& destination_id_map,
                  handle_t source_handle){

    auto id = source_graph.get_id(source_handle);
    auto name = source_id_map.get_name(id);

    nid_t translated_id;
    if (destination_id_map.exists(name)){
        translated_id = destination_id_map.get_id(name);
    }
    else{
        translated_id = destination_id_map.insert(name);
    }

    return {id,translated_id};
}


void split_connected_components(
        MutablePathDeletableHandleGraph& graph,
        IncrementalIdMap<string>& id_map,
        vector<HashGraph>& graphs,
        vector<IncrementalIdMap<string> >& id_maps,
        bool delete_visited_components) {

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
        unordered_set<nid_t> to_be_deleted;

        auto start_node = *all_nodes.begin();

        // Duplicate all the nodes
        for_node_in_bfs(graph, start_node, [&](const handle_t& h) {
            auto s = graph.get_sequence(h);
            nid_t id;
            nid_t other_id;
            tie(id, other_id) = translate_id(graph, id_map, graphs.back(), id_maps.back(), h);

            cerr << "iterating " << id << '\n';

            assert(not graphs.back().has_node(other_id));
            graphs.back().create_handle(s, other_id);

            graph.for_each_step_on_handle(h, [&](const step_handle_t s){
                paths_to_be_copied.emplace(graph.get_path_name(graph.get_path_handle_of_step(s)));
            });

            if (delete_visited_components) {
                to_be_deleted.emplace(id);
            }

            if (id != start_node) {
                all_nodes.erase(id);
                cerr << "erasing " << id << '\n';
            }
        });

        // Duplicate all the edges
        for_edge_in_bfs(graph, start_node, [&](const handle_t& handle_a, const handle_t& handle_b) {
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

        all_nodes.erase(start_node);
        cerr << "erasing " << start_node << '\n';

        if (delete_visited_components) {
            for (auto& n: to_be_deleted) {
                cerr << "deleting " << n << '\n';
                auto h = graph.get_handle(n);
                graph.destroy_handle(h);
            }
        }
    }
}


void write_connected_components_to_gfas(
        const MutablePathDeletableHandleGraph& graph,
        const IncrementalIdMap<string>& id_map,
        path output_directory) {

    unordered_set<nid_t> all_nodes;

    graph.for_each_handle([&](const handle_t& h) {
        auto n = graph.get_id(h);
        all_nodes.emplace(n);
    });

    size_t i = 0;
    while (not all_nodes.empty()) {
        cerr << "NEW CONNECTED COMPONENT" << '\n';
        string filename_prefix = output_directory / ("component_" + to_string(i));
        ofstream file(filename_prefix + ".gfa");

        unordered_set<string> paths_to_be_copied;
        unordered_set<nid_t> to_be_deleted;

        auto start_node = *all_nodes.begin();

        // Duplicate all the nodes
        for_node_in_bfs(graph, start_node, [&](const handle_t& h) {
            auto id = graph.get_id(h);
            if (id != start_node){
                all_nodes.erase(id);
            }

            write_node_to_gfa(graph, id_map, h, file);

            graph.for_each_step_on_handle(h, [&](const step_handle_t s){
                paths_to_be_copied.emplace(graph.get_path_name(graph.get_path_handle_of_step(s)));
            });

        });

        // Duplicate all the edges
        for_edge_in_bfs(graph, start_node, [&](const handle_t& handle_a, const handle_t& handle_b) {
            write_edge_to_gfa(graph, id_map, {handle_a, handle_b}, file);
        });

        // Duplicate all the paths
        for (auto& path_name: paths_to_be_copied){
            write_path_to_gfa(graph, id_map, graph.get_path_handle(path_name), file);
        }

        all_nodes.erase(start_node);
        i++;
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

// Find any nodes that are adjacent to the beginning and end of a path, as long as they are the only adjacent node
pair<handle_t, bool> find_singleton_adjacent_handle(const PathHandleGraph& graph, const handle_t& h, bool left) {
    handle_t adjacent_handle;
    size_t n_adjacent = 0;
    bool success = false;

    graph.follow_edges(h, left, [&](const handle_t& other_handle) {
        if (n_adjacent > 0) {
            return false;   // Exit loop
        }

        adjacent_handle = other_handle;
        n_adjacent++;

        return true;
    });

    // Only return true if there was exactly one adjacent handle
    if (n_adjacent == 1) {
        success = true;
    }

    return {adjacent_handle, success};
}


/// Cheap way to check if path is part of a diploid phased pair of paths. It actually just relies on
/// the Shasta convention that phased paths always end on a bubble
void find_diploid_paths(const PathHandleGraph& graph, vector<path_handle_t>& diploid_paths){
    graph.for_each_path_handle([&](const path_handle_t& p){
        auto begin_handle = graph.get_handle_of_step(graph.path_begin(p));
        auto end_handle = graph.get_handle_of_step(graph.path_back(p));

        bool begin_is_bubble = find_singleton_adjacent_handle(graph, begin_handle, true).second;
        bool end_is_bubble = find_singleton_adjacent_handle(graph, end_handle, false).second;

        if (begin_is_bubble and end_is_bubble){
            diploid_paths.emplace_back(p);
        }
    });
}


/// Cheap way to check if path is part of a diploid phased pair of paths. It actually just relies on
/// the Shasta convention that phased paths always end on a bubble
void find_diploid_paths(const PathHandleGraph& graph, const set<string>& subset, vector<path_handle_t>& diploid_paths, char path_delimiter){
    vector<path_handle_t> paths;

    graph.for_each_path_handle([&](const path_handle_t& p){
        string name = graph.get_path_name(p);

        string graph_component;
        size_t component_haplotype;

        tie(graph_component, component_haplotype) = parse_path_string(name, path_delimiter);

        if (subset.find(graph_component) != subset.end()){
            paths.emplace_back(p);
        }
    });

    for (const auto& p: paths){
        auto begin_handle = graph.get_handle_of_step(graph.path_begin(p));
        auto end_handle = graph.get_handle_of_step(graph.path_back(p));

        bool begin_is_bubble = find_singleton_adjacent_handle(graph, begin_handle, true).second;
        bool end_is_bubble = find_singleton_adjacent_handle(graph, end_handle, false).second;

        if (begin_is_bubble and end_is_bubble){
            diploid_paths.emplace_back(p);
        }
    }
}


void extend_paths(
        MutablePathMutableHandleGraph& graph,
        vector<pair<path_handle_t, handle_t> >& to_be_prepended,
        vector<pair<path_handle_t, handle_t> >& to_be_appended) {

    graph.for_each_path_handle([&](const path_handle_t& p) {
        auto begin_handle = graph.get_handle_of_step(graph.path_begin(p));
        auto end_handle = graph.get_handle_of_step(graph.path_back(p));

        auto left_result = find_singleton_adjacent_handle(graph, begin_handle, true);
        auto right_result = find_singleton_adjacent_handle(graph, end_handle, false);

        if (left_result.second) {
            to_be_prepended.emplace_back(p, left_result.first);
        }

        if (right_result.second) {
            to_be_appended.emplace_back(p, right_result.first);
        }
    });

    for (auto& item: to_be_prepended) {
        // PREpend the LEFT side node if it meets the conditions
        graph.prepend_step(item.first, item.second);
    }

    for (auto& item: to_be_appended) {
        // Append the RIGHT side node if it meets the conditions
        graph.append_step(item.first, item.second);
    }
}


void un_extend_paths(
        MutablePathMutableHandleGraph& graph,
        const vector <pair<path_handle_t, handle_t> >& to_be_prepended,
        const vector <pair<path_handle_t, handle_t> >& to_be_appended) {

    for (auto& item: to_be_prepended) {
        // un-PREpend (pop) the LEFT side node
        auto begin = graph.path_begin(item.first);
        auto next = graph.get_next_step(begin);
        graph.rewrite_segment(begin, next, vector<handle_t>());
    }

    for (auto& item: to_be_appended) {
        // un-apppend (pop) the RIGHT side node
        auto end = graph.path_end(item.first);
        auto prev = graph.get_previous_step(end);
        graph.rewrite_segment(prev, end, vector<handle_t>());
    }
}


void unzip(MutablePathDeletableHandleGraph& graph, IncrementalIdMap<string>& id_map){
    unordered_set<nid_t> nodes_to_be_destroyed;

    cout << graph.get_path_count() << '\n';
    vector<string> path_names;

    vector<path_handle_t> paths;

    graph.for_each_path_handle([&](const path_handle_t& p) {
        paths.emplace_back(p);
    });

    for (auto& p: paths){
        string path_sequence;
        cerr << "Path " << graph.get_path_name(p) << '\n';

        graph.for_each_step_in_path(p, [&](const step_handle_t s){
            handle_t h = graph.get_handle_of_step(s);
            nid_t n = graph.get_id(h);

            string sequence = graph.get_sequence(h);
            path_sequence += sequence;

            string name = id_map.get_name(n);
            cerr << '\t' << name << " " << (sequence.size() < 100 ? sequence : "--") << '\n';

            nodes_to_be_destroyed.emplace(n);
        });

        string haplotype_path_name = graph.get_path_name(p) + "_hap";
        auto new_id = id_map.insert(haplotype_path_name);
        handle_t haplotype_handle = graph.create_handle(path_sequence, new_id);

        auto path_start_handle = graph.get_handle_of_step(graph.path_begin(p));
        auto path_stop_handle = graph.get_handle_of_step(graph.path_back(p));

        // Find neighboring nodes for the path and create edges to the new haplotype node (LEFT)
        graph.follow_edges(path_start_handle, true, [&](const handle_t& other){
            graph.create_edge(other, haplotype_handle);
        });

        // Find neighboring nodes for the path and create edges to the new haplotype node (RIGHT)
        graph.follow_edges(path_stop_handle, false, [&](const handle_t& other){
            graph.create_edge(haplotype_handle, other);
        });

        // Label the new haplotype node using a path_name to indicate which path it is derived from
        // TODO: track provenance, update id_map?
        auto haplotype_path_handle = graph.create_path_handle(haplotype_path_name);
        graph.append_step(haplotype_path_handle, haplotype_handle);
    }

    // Destroy the nodes that have had their sequences duplicated into haplotypes
    for (auto& n: nodes_to_be_destroyed){
        cerr << "Destroying: " << id_map.get_name(n) << '\n';
        graph.destroy_handle(graph.get_handle(n));
    }

    // TODO: when creating haplotypes, maintain error bubbles?
    // Search for islands that were created by unphased nodes when the haplotype paths were deleted,
    // and delete the islands (assuming they were errors, and downstream applications won't want bubbles)
    for_each_connected_component(graph, [&](unordered_set<nid_t>& component){
        bool has_path = false;

        // Iterate the connected component and check if it has any haplotype info (if not, it was a unlabeled bubble)
        for (auto& id: component){
            auto h = graph.get_handle(id);

            graph.for_each_step_on_handle(h, [&](const step_handle_t s){
                has_path = true;
                return false;
            });

//            cerr << id_map.get_name(graph.get_id(h)) << " has_path: " << has_path << '\n';
        }

        // Delete any component that has no path label
        if (not has_path) {
            for (auto& id: component) {
                auto h = graph.get_handle(id);
                graph.destroy_handle(h);
            }
        }
    });
}



}

