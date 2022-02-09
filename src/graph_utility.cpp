#include "graph_utility.hpp"

namespace gfase {


pair<string, size_t> parse_path_string(string path_name, char delimiter){
    size_t index = path_name.rfind(delimiter);
    string component_name = path_name.substr(0,index);
    size_t component_haplotype = stoi(path_name.substr(index+1,path_name.length()));

    return {component_name, component_haplotype};
}


void for_node_in_bfs(
        const HandleGraph& graph,
        nid_t start_node,
        const function<bool(const nid_t& id)>& pass_criteria,
        const function<void(const handle_t& h)>& f){

    // Do nothing if the start node was blacklisted
    if (not pass_criteria(start_node)){
        return;
    }

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
            auto pass_a = visited.emplace(other_node).second;

            // Check if this node is excluded by choice of the user
            auto pass_b = pass_criteria(other_node);

            // Check that this has NOT been visited before queuing it
            if (pass_a and pass_b) {
                q.emplace(other_node);
            }
        });

        graph.follow_edges(h, true, [&](const handle_t& other_handle) {
            auto other_node = graph.get_id(other_handle);

            // Attempt to add this to the set of visited nodes
            auto pass_a = visited.emplace(other_node).second;

            // Check if this node is excluded by choice of the user
            auto pass_b = pass_criteria(other_node);

            // Check that this has NOT been visited before queuing it
            if (pass_a and pass_b) {
                q.emplace(other_node);
            }
        });
    }
}


void for_node_in_bfs(
        const HandleGraph& graph,
        nid_t start_node,
        const unordered_set<nid_t>& do_not_visit,
        const function<void(const handle_t&)>& f){

    for_node_in_bfs(
            graph,
            start_node,
            [&](const nid_t& id){return (do_not_visit.find(id) == do_not_visit.end());},
            f);
}

void for_node_in_bfs(const HandleGraph& graph, nid_t start_node, const function<void(const handle_t&)>& f) {
    unordered_set<nid_t> do_not_visit;
    for_node_in_bfs(graph, start_node, do_not_visit, f);
}


void for_edge_in_bfs(
        const HandleGraph& graph,
        nid_t start_node,
        const unordered_set<nid_t>& do_not_visit,
        const function<void(const handle_t& handle_a, const handle_t& handle_b)>& f_pass,
        const function<void(const handle_t& handle_a, const handle_t& handle_b)>& f_fail) {

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
            auto pass_a = visited_nodes.emplace(other_node).second;

            // Check if this node has been blacklisted
            auto pass_b = do_not_visit.find(other_node) == do_not_visit.end();

            // Check that this has NOT been visited before queuing it
            if (pass_a and pass_b) {
                q.emplace(other_handle);
            }

            auto edge_result = visited_edges[h].insert(other_handle);

            if (edge_result.second){
                if (pass_b) {
                    f_pass(h, other_handle);
                }
                else {
                    f_fail(h, other_handle);
                }
            }
        });

        graph.follow_edges(h, true, [&](const handle_t& other_handle) {
            auto other_node = graph.get_id(other_handle);

            // Attempt to add this to the set of visited nodes
            auto pass_a = visited_nodes.emplace(other_node).second;

            // Check if this node has been blacklisted
            auto pass_b = do_not_visit.find(other_node) == do_not_visit.end();

            // Check that this has NOT been visited before queuing it
            if (pass_a and pass_b) {
                q.emplace(other_handle);
            }

            auto edge_result = visited_edges[other_handle].emplace(h);

            if (edge_result.second){
                if (pass_b) {
                    f_pass(other_handle, h);
                }
                else {
                    f_fail(other_handle, h);
                }
            }
        });

        if (begin){
            begin = false;
        }
    }
}


// Simplifying wrapper for the "for_edge_in_bfs" method, which has some extra criteria which is not always necessary
void for_edge_in_bfs(
        const HandleGraph& graph,
        nid_t start_node,
        const function<void(const handle_t& handle_a, const handle_t& handle_b)>& f) {

    unordered_set<nid_t> do_not_visit;
    for_edge_in_bfs(graph, start_node, do_not_visit, f, [&](const handle_t& handle_a, const handle_t& handle_b){});
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


tuple<nid_t,nid_t,bool> try_translate_id(const HandleGraph& source_graph,
                               const IncrementalIdMap<string>& source_id_map,
                               HandleGraph& destination_graph,
                               IncrementalIdMap<string>& destination_id_map,
                               handle_t source_handle){

    auto id = source_graph.get_id(source_handle);
    auto name = source_id_map.get_name(id);
    bool success = destination_id_map.exists(name);

    nid_t translated_id;
    if (success){
        translated_id = destination_id_map.get_id(name);
    }

    return {id, translated_id, success};
}


void split_connected_components(
        MutablePathDeletableHandleGraph& graph,
        IncrementalIdMap<string>& id_map,
        vector<HashGraph>& graphs,
        vector<IncrementalIdMap<string> >& id_maps,
        vector <vector <pair <string, string> > >& in_edges,
        vector <vector <pair <string, string> > >& out_edges,
        const unordered_set<nid_t>& do_not_visit,
        bool delete_visited_components) {

    unordered_set<nid_t> all_nodes;

    graph.for_each_handle([&](const handle_t& h) {
        auto n = graph.get_id(h);
        all_nodes.emplace(n);
    });

    while (not all_nodes.empty()) {
        auto start_node = *all_nodes.begin();

//        cerr << "Trying: " << id_map.get_name(start_node) << '\n';

        // Skip the node if it is blacklisted
        if (do_not_visit.find(start_node) != do_not_visit.end()){
            all_nodes.erase(start_node);
            continue;
        }

        // Allocate new elements in the vectors for this component
        graphs.emplace_back();
        id_maps.emplace_back();
        in_edges.emplace_back();
        out_edges.emplace_back();

        unordered_set<string> paths_to_be_copied;
        unordered_set<nid_t> to_be_deleted;

        // Duplicate all the nodes
        for_node_in_bfs(graph, start_node, do_not_visit, [&](const handle_t& h) {
            auto s = graph.get_sequence(h);
            nid_t id;
            nid_t other_id;
            tie(id, other_id) = translate_id(graph, id_map, graphs.back(), id_maps.back(), h);

//            cerr << "Iterating: " << id_map.get_name(id) << '\n';

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
            }
        });

        // Duplicate all the edges
        for_edge_in_bfs(graph, start_node, do_not_visit, [&](const handle_t& handle_a, const handle_t& handle_b) {
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
        },
        [&](const handle_t& handle_a, const handle_t& handle_b){
            auto id_a = graph.get_id(handle_a);
            auto id_b = graph.get_id(handle_b);
            auto name_a = id_map.get_name(id_a);
            auto name_b = id_map.get_name(id_b);

            if (do_not_visit.find(id_a) != do_not_visit.end()){
                in_edges.back().emplace_back(name_a, name_b);
            }
            else if (do_not_visit.find(id_b) != do_not_visit.end()){
                out_edges.back().emplace_back(name_a, name_b);
            }
            else{
                throw runtime_error("ERROR: unexpected membership in do_not_visit for edge: " + name_a + " -> " + name_b);
            }
        });

        // Duplicate all the paths
        for (auto& path_name: paths_to_be_copied){
            auto p = graph.get_path_handle(path_name);

            assert(not graphs.back().has_path(path_name));
            path_handle_t other_p;
            bool prev_success = false;
            size_t n_breaks = 0;

            graph.for_each_step_in_path(p, [&](const step_handle_t& s){
                auto h = graph.get_handle_of_step(s);

                // Check if this node exists in the new subgraph, and if so, find its ID
                nid_t id;
                nid_t other_id;
                bool success;
                tie(id, other_id, success) = try_translate_id(graph, id_map, graphs.back(), id_maps.back(), h);

                // It is possible that the path has been broken during BFS, if a node has been (optionally) excluded
                // by the user. In that case, make a new path for every time there is a break in path continuity.
                if (success) {
                    if (not prev_success) {
                        string new_path_name = path_name + (n_breaks > 0 ? to_string(n_breaks) : "");
                        other_p = graphs.back().create_path_handle(path_name);
                        n_breaks++;
                    }

                    auto other_h = graph.get_handle(other_id, graph.get_is_reverse(h));
                    graphs.back().append_step(other_p, other_h);
                }

                prev_success = success;
            });
        }

        all_nodes.erase(start_node);

        if (delete_visited_components) {
            for (auto& n: to_be_deleted) {
                auto h = graph.get_handle(n);
                graph.destroy_handle(h);
            }
        }
    }
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

        if (graphs.back().get_node_count() == 0){
            throw runtime_error("connected component with start name " + id_map.get_name(start_node) + " creates empty graph");
        }

        all_nodes.erase(start_node);

        if (delete_visited_components) {
            for (auto& n: to_be_deleted) {
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
void find_diploid_paths(const PathHandleGraph& graph, unordered_set<string>& diploid_path_names){
    graph.for_each_path_handle([&](const path_handle_t& p){
        auto begin_handle = graph.get_handle_of_step(graph.path_begin(p));
        auto end_handle = graph.get_handle_of_step(graph.path_back(p));

        bool begin_is_bubble = find_singleton_adjacent_handle(graph, begin_handle, true).second;
        bool end_is_bubble = find_singleton_adjacent_handle(graph, end_handle, false).second;

        if (begin_is_bubble and end_is_bubble){
            diploid_path_names.emplace(graph.get_path_name(p));
        }
    });
}


/// Exhaustive check for overlapping paths, building a bidirectional mapping of diploid paths to one another
void find_diploid_paths(
        const PathHandleGraph& graph,
        unordered_map<string, string>& diploid_path_names,
        unordered_set<string>& haploid_path_names){

    graph.for_each_path_handle([&](const path_handle_t& p){
        auto path_name = graph.get_path_name(p);

        unordered_set <path_handle_t> overlapping_paths;
        graph.for_each_step_in_path(p, [&](const step_handle_t& s){
            auto h = graph.get_handle_of_step(s);

            graph.for_each_step_on_handle(h, [&](const step_handle_t& s_other){
                if (graph.get_path_handle_of_step(s_other) != p) {
                    overlapping_paths.emplace(graph.get_path_handle_of_step(s_other));
                }
            });
        });

        if (overlapping_paths.empty()){
            haploid_path_names.emplace(path_name);
        }
        else if (overlapping_paths.size() == 1){
            string other_path_name = graph.get_path_name(*overlapping_paths.begin());

            // Build mapping in both directions
            diploid_path_names[other_path_name] = path_name;
            diploid_path_names[path_name] = other_path_name;
        }
        else{
            cerr << "Found more than 1 overlapping path\n\t";
            for (auto& p_other: overlapping_paths){
                cerr << graph.get_path_name(p_other) << ' ';
            }
            cerr << '\n';

            throw runtime_error("ERROR: path overlaps with more than one other (is not diploid): " + path_name);
        }

    });
}


/// Cheap way to check if path is part of a diploid phased pair of paths. It actually just relies on
/// the Shasta convention that phased paths always end on a bubble
void find_diploid_paths(const PathHandleGraph& graph, const set<string>& subset, unordered_set<string>& diploid_path_names, char path_delimiter){
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
            diploid_path_names.emplace(graph.get_path_name(p));
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
        // un-append (pop) the RIGHT side node
        auto end = graph.path_end(item.first);
        auto prev = graph.get_previous_step(end);
        graph.rewrite_segment(prev, end, vector<handle_t>());
    }
}


void unzip(MutablePathDeletableHandleGraph& graph, IncrementalIdMap<string>& id_map, bool keep_paths){
    unordered_set<nid_t> nodes_to_be_destroyed;
    vector<string> path_names;
    string temporary_suffix = "_hap";

    vector<path_handle_t> paths;

//    cerr << "Paths in component:" << '\n';
    graph.for_each_path_handle([&](const path_handle_t& p) {
        paths.emplace_back(p);
        path_names.emplace_back(graph.get_path_name(p));
//        cerr << graph.get_path_name(p) << '\n';
    });

    for (auto& p: paths){
        string path_sequence;

        graph.for_each_step_in_path(p, [&](const step_handle_t s){
            handle_t h = graph.get_handle_of_step(s);
            nid_t n = graph.get_id(h);

            string sequence = graph.get_sequence(h);
            path_sequence += sequence;

            string name = id_map.get_name(n);

            nodes_to_be_destroyed.emplace(n);
        });

        auto name = graph.get_path_name(p);
        string haplotype_path_name = name + temporary_suffix;

        cerr << haplotype_path_name << '\n';

        auto new_id = id_map.insert(name);
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
        graph.destroy_handle(graph.get_handle(n));
    }

    // If there were no paths in this entire component, dont delete anything, just leave it as is
    if (paths.empty()){
        return;
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

        }

        // Delete any component that has no path label
        if (not has_path) {
            for (auto& id: component) {
                auto h = graph.get_handle(id);
                graph.destroy_handle(h);
            }
        }
    });

    // Go back and rename all the paths so they don't have the unnecessary added suffix
    for (const auto& name: path_names){
        // Find the path that was renamed with a suffix, and find its only node (all unzipped paths are singletons)
        auto path_name_with_suffix = name + temporary_suffix;
        auto p_suffix = graph.get_path_handle(path_name_with_suffix);

        if (keep_paths) {
            auto h = graph.get_handle_of_step(graph.path_begin(p_suffix));

            // Make a copy without the suffix
            auto p = graph.create_path_handle(name);
            graph.append_step(p, h);
        }

        // Destroy the path with the suffix
        graph.destroy_path(p_suffix);
    }
}


void for_each_tip(const HandleGraph& graph, const function<void(const handle_t& h, bool is_left, bool is_right)>& f){
    graph.for_each_handle([&](const handle_t& h_i){
        bool is_left = graph.get_degree(h_i,true);
        bool is_right = graph.get_degree(h_i,false);

        f(h_i, is_left, is_right);
    });
}


void write_paths_to_csv(const PathHandleGraph& graph, const IncrementalIdMap<string>& id_map, ofstream& file){
    graph.for_each_path_handle([&](const path_handle_t& path){
        string path_name = graph.get_path_name(path);
        size_t n_steps = graph.get_step_count(path);
        size_t i = 0;

        file << path_name << ',' << n_steps << ',';

        graph.for_each_step_in_path(path, [&](const step_handle_t& s){
            auto h = graph.get_handle_of_step(s);
            auto name = id_map.get_name(graph.get_id(h));

            file << name << (graph.get_is_reverse(h) ? '-' : '+');
            if (i < n_steps - 1){
                file << ' ';
            }

            i++;
        });

        file << '\n';
    });
}


}

