#include "Bridges.hpp"

#include <functional>
#include <unordered_map>
#include <utility>

using std::function;
using std::pair;

using handlegraph::as_integer;

namespace gfase {

/*
 * Wrapper for a handle graph that uses the bi-edged formalism instead
 * of a bidirected formalism, which allows a simple undirected graph
 * interface
 */
class BiedgedGraph {
public:
    BiedgedGraph(const HandleGraph& graph);
    BiedgedGraph() = delete;
    
    size_t size() const;
    
    // returns an invalid handle if graph is empty
    handle_t get_arbitrary_node() const;
    
    void for_each_node(const function<void(handle_t)>& lambda) const;
    // bool is true for across-node edges and false for adjacency edges
    void for_each_neighbor(handle_t node, const function<void(handle_t, bool)>& lambda) const;
    
private:
    
    const HandleGraph& graph;
    
};

BiedgedGraph::BiedgedGraph(const HandleGraph& graph) : graph(graph) {
    
}

size_t BiedgedGraph::size() const {
    return 2 * graph.size();
}

handle_t BiedgedGraph::get_arbitrary_node() const {
    handle_t handle;
    as_integer(handle) = -1; // to avoid uninitialized warnings
    graph.for_each_handle([&handle](const handle_t& h) {
        handle = h;
        return false;
    });
    return handle;
}


void BiedgedGraph::for_each_node(const function<void(handle_t)>& lambda) const {
    graph.for_each_handle([&graph](const handle_t& h) {
        lambda(h);
        lambda(graph.flip(h));
    });
}

void BiedgedGraph::for_each_neighbor(handle_t node, const function<void(handle_t, bool)>& lambda) const {
    
    lambda(graph.flip(node), true);
    graph.follow_edges(node, false, [&graph, &lambda](const handle_t& next) {
        lambda(graph.flip(next), false);
    });
}

vector<handle_t> bridge_nodes(const HandleGraph& graph) {
    
    BiedgedGraph undir_graph(graph);
    
    // TODO: should I just make this a class?
    // records of node -> (edge upward, traversed upward, dfs index, edges downward)
    unordered_map<handle_t, tuple<pair<handle_t, bool>, bool, size_t, vector<pair<handle_t, bool>>>> dfs_tree;
    vector<handle_t> dfs_order;
    
    dfs_tree.reserve(undir_graph.size());
    dfs_order.reserve(undir_graph.size());
    
    // records of (node, next edge idx, edges)
    vector<tuple<handle_t, size_t, vector<pair<handle_t, bool>>>> stack;
    
    // enqueue a non-visited node
    auto enqueue = [&](handle_t node, handle_t parent, bool from_across_edge) {
        
        stack.emplace_back(node, 0, vector<pair<handle_t, bool>>());
        
        // record the order and initialize the tree edges
        dfs_tree[node].emplace_back({parent, from_across_edge}, false, dfs_order.sizse(), {});
        dfs_order.emplace_back(node);
        
        undir_graph.for_each_neighbor(node, [&](handle_t nbr, bool is_across_edge) {
            get<2>(stack.back()).emplace_back(nbr, is_across_edge);;
        });
    };
    
    // DFS
    handle_t root = undir_graph.get_arbitrary_node();
    enqueue(root, root, false);
    // note: the root doesn't really have a parent, but we use itself as a sentinel
    while (!stack.empty()) {
        
        auto& top = stack.back();
        
        if (get<1>(top) == get<2>(top).size()) {
            // exhausted edges from this node
            stack.pop_back();
            continue;
        }
        
        auto next = get<2>(top)[get<1>(top)++];
        
        if (!dfs_tree.count(next.first)) {
            // the next node hasn't been visited yet, so we add its edge to
            // the DFS tree and queue it up
            get<3>(dfs_tree[get<0>(top)]).emplace_back(next);
            enqueue(next.first, get<0>(top), next.second);
        }
    }
    
    for (size_t i = 0; i < dfs_order.size(); ++i) {
        
        handle_t node = dfs_order[i];
        
        // get the non-tree edges that will initiate chains
        vector<handle_t> chain_heads;
        // note: we count on the edges being iterated in the same order as previously
        const auto& tree_record = dfs_order[node];
        size_t j = 0;
        undir_graph.for_each_neighbor(node, [&](handle_t nbr, bool is_across_edge) {
            if (j < get<3>(tree_record).size() && get<3>(tree_record)[j] == make_pair(nbr, is_across_edge)) {
                // this edge was used in the DFS tree
                ++j;
            }
            else if (get<2>(tree_record) > i) {
                // this edge points down the tree rather than up it
                chain_heads.emplace_back(nbr);
            }
        });
        
        // mark nodes in the untraversed portion of this chain's cycle as traversed
        for (handle_t chain_cursor : chain_heads) {
            while (chain_cursor != node && !get<1>(dfs_tree[chain_cursor])) {
                get<1>(dfs_tree[chain_cursor]) = true;
                chain_cursor = get<0>(dfs_tree[chain_cursor]).first;
            }
        }
    }
    
    // find the across-node edges that were never traversed in a cycle
    vector<handle_t> bridges;
    
    for (const auto& tree_record : dfs_tree) {
        if (!get<1>(tree_record.second) && get<0>(tree_record.second).second) {
            bridges.push_back(graph.forward(tree_record.first));
        }
    }
    
    return bridges;
}

}
