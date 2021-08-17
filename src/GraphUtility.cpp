#include "GraphUtility.hpp"

namespace gfase{

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
            if (result.second){
                q.emplace(other_node);
            }
        });

        graph.follow_edges(h, true, [&](const handle_t& other_handle) {
            auto other_node = graph.get_id(other_handle);

            // Attempt to add this to the set of visited nodes
            auto result = visited.emplace(other_node);

            // Check that this has NOT been visited before queuing it
            if (result.second){
                q.emplace(other_node);
            }
        });
    }
}

}