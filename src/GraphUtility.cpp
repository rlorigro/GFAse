#include "GraphUtility.hpp"

namespace gfase {

void for_node_in_bfs(HandleGraph& graph, nid_t start_node, const function<void(const handle_t&)>& f) {
    unordered_set<nid_t> visited;
    queue<nid_t> q;

    q.emplace(start_node);

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
            visited_nodes.emplace(start_node);
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


pair<nid_t,nid_t> translate_id(HandleGraph& source_graph,
                  IncrementalIdMap<string>& source_id_map,
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


void split_connected_components(HandleGraph& graph, IncrementalIdMap<string>& id_map, vector<HashGraph>& graphs, vector<IncrementalIdMap<string> >& id_maps) {
    unordered_set<nid_t> all_nodes;

    graph.for_each_handle([&](const handle_t& h) {
        auto n = graph.get_id(h);
        all_nodes.emplace(n);
    });

    while (not all_nodes.empty()) {
        graphs.emplace_back();
        id_maps.emplace_back();

        auto iter = all_nodes.begin();

        for_node_in_bfs(graph, *iter, [&](const handle_t& h) {
            auto s = graph.get_sequence(h);
            nid_t id;
            nid_t other_id;
            tie(id, other_id) = translate_id(graph, id_map, graphs.back(), id_maps.back(), h);

            auto other_handle = graphs.back().create_handle(s, other_id);

            all_nodes.erase(id);
        });

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

            graphs.back().create_edge(other_handle_a, other_handle_b);
        });


    }
}


}

