#include "Bridges.hpp"

#include "SubgraphOverlay.hpp"

#include "handlegraph/util.hpp"

#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <utility>
#include <iostream>

using std::tuple;
using std::get;
using std::tie;
using std::make_pair;
using std::make_tuple;
using std::reverse;
using std::unordered_map;
using std::unordered_set;
using std::cerr;
using std::endl;

using handlegraph::as_integer;
using handlegraph::as_handle;
using handlegraph::nid_t;

namespace gfase {

static const bool debug = false;

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
    void for_each_neighbor(handle_t node, const function<void(handle_t)>& lambda) const;
    bool is_across_node_edge(handle_t node1, handle_t node2) const;
    
private:
    
    const HandleGraph& graph;
    
};

BiedgedGraph::BiedgedGraph(const HandleGraph& graph) : graph(graph) {
    
}

size_t BiedgedGraph::size() const {
    return 2 * graph.get_node_count();
}

handle_t BiedgedGraph::get_arbitrary_node() const {
    handle_t handle = as_handle(-1); // to avoid uninitialized warnings
    graph.for_each_handle([&handle](const handle_t& h) {
        handle = h;
        return false;
    });
    return handle;
}

void BiedgedGraph::for_each_node(const function<void(handle_t)>& lambda) const {
    graph.for_each_handle([&](const handle_t& h) {
        lambda(h);
        lambda(graph.flip(h));
    });
}

void BiedgedGraph::for_each_neighbor(handle_t node, const function<void(handle_t)>& lambda) const {
    
    lambda(graph.flip(node));
    graph.follow_edges(node, false, [&](const handle_t& next) {
        lambda(graph.flip(next));
    });
}

bool BiedgedGraph::is_across_node_edge(handle_t node1, handle_t node2) const {
    return graph.flip(node1) == node2;
}

vector<handle_t> bridge_nodes(const HandleGraph& graph) {
    
    if (debug) {
        cerr << "finding bridges in graph" << endl;
        graph.for_each_handle([&](const handle_t& handle) {
            cerr << graph.get_id(handle) << " " << graph.get_sequence(handle) << endl;
            graph.follow_edges(handle, true, [&](const handle_t& prev) {
                cerr << "\t" << graph.get_id(prev) << (graph.get_is_reverse(prev) ? "-" : "+") << " <-" << endl;
            });
            graph.follow_edges(handle, false, [&](const handle_t& next) {
                cerr << "\t-> " << graph.get_id(next)  << (graph.get_is_reverse(next) ? "-" : "+")<< endl;
            });
        });
    }
    
    BiedgedGraph undir_graph(graph);
    
    // TODO: should I just make this a class?
    // records of node -> (parent, traversed upward, dfs index, edges downward)
    unordered_map<handle_t, tuple<handle_t, bool, size_t, vector<handle_t>>> dfs_tree;
    vector<handle_t> dfs_order;
    
    dfs_tree.reserve(undir_graph.size());
    dfs_order.reserve(undir_graph.size());
    
    // records of (node, next edge idx, edges)
    vector<tuple<handle_t, size_t, vector<handle_t>>> stack;
    
    // enqueue a non-visited node
    auto enqueue = [&](handle_t node, handle_t parent) {
        
        stack.emplace_back(node, 0, vector<handle_t>());
        
        // record the order and initialize the tree edges
        dfs_tree[node] = make_tuple(parent, false, dfs_order.size(), vector<handle_t>());
        dfs_order.emplace_back(node);
        
        undir_graph.for_each_neighbor(node, [&](handle_t nbr) {
            get<2>(stack.back()).emplace_back(nbr);
        });
    };
    
    // DFS
    handle_t root = undir_graph.get_arbitrary_node();
    if (debug) {
        cerr << "beginning DFS from " << graph.get_id(root) << (graph.get_is_reverse(root) ? "L" : "R") << endl;
    }
    enqueue(root, root);
    // note: the root doesn't really have a parent, but we use itself as a sentinel
    while (!stack.empty()) {
        
        auto& top = stack.back();
        
        if (get<1>(top) == get<2>(top).size()) {
            // exhausted edges from this node
            stack.pop_back();
            continue;
        }
        
        auto next = get<2>(top)[get<1>(top)++];
        
        if (!dfs_tree.count(next)) {
            // the next node hasn't been visited yet, so we add its edge to
            // the DFS tree and queue it up
            get<3>(dfs_tree[get<0>(top)]).emplace_back(next);
            enqueue(next, get<0>(top));
        }
    }
    
    if (debug) {
        cerr << "finished making DFS tree:" << endl;
        for (const auto& node : dfs_tree) {
            cerr << graph.get_id(node.first) << (graph.get_is_reverse(node.first) ? "L" : "R") << endl;
            for (auto next : get<3>(node.second)) {
                cerr << "\t-> " << graph.get_id(next) << (graph.get_is_reverse(next) ? "L" : "R") << endl;
            }
        }
    }
    
    for (size_t i = 0; i < dfs_order.size(); ++i) {
        
        handle_t node = dfs_order[i];
        
        // get the non-tree edges that will initiate chains
        vector<handle_t> chain_heads;
        // note: we count on the edges being iterated in the same order as previously
        const auto& tree_record = dfs_tree[node];
        size_t j = 0;
        if (debug) {
            cerr << "checking for unused edges from " << graph.get_id(node) << (graph.get_is_reverse(node) ? "L" : "R") << endl;
        }
        undir_graph.for_each_neighbor(node, [&](handle_t nbr) {
            if (debug) {
                cerr << "\t" << graph.get_id(nbr) << (graph.get_is_reverse(nbr) ? "L" : "R");
            }
            if (j < get<3>(tree_record).size() && get<3>(tree_record)[j] == nbr) {
                // this edge was used in the DFS tree
                ++j;
                if (debug) {
                    cerr << " was used in tree" << endl;
                }
            }
            else if (get<2>(dfs_tree[nbr]) > i) {
                // this edge points down the tree rather than up it
                chain_heads.emplace_back(nbr);
                if (debug) {
                    cerr << " is a downward edge" << endl;
                }
            }
            else {
                if (debug) {
                    cerr << " is an upward edge" << endl;
                }
            }
        });
        
        // mark nodes in the untraversed portion of this chain's cycle as traversed
        for (handle_t chain_cursor : chain_heads) {
            if (debug) {
                cerr << "beginning upward traversal from " << graph.get_id(chain_cursor) << (graph.get_is_reverse(chain_cursor) ? "L" : "R") << endl;
            }
            while (chain_cursor != node && !get<1>(dfs_tree[chain_cursor])) {
                get<1>(dfs_tree[chain_cursor]) = true;
                chain_cursor = get<0>(dfs_tree[chain_cursor]);
            }
        }
    }
    
    // find the across-node edges that were never traversed in a cycle
    vector<handle_t> bridges;
    
    for (const auto& tree_record : dfs_tree) {

        if (!get<1>(tree_record.second) &&
            undir_graph.is_across_node_edge(tree_record.first, get<0>(tree_record.second))) {
            if (debug) {
                cerr << "marking edge " << graph.get_id(tree_record.first) << " - " << (graph.get_is_reverse(tree_record.first) ? "L" : "R") << " as a bridge node" << endl;
            }
            bridges.push_back(graph.forward(tree_record.first));
        }
    }
    
    return bridges;
}

vector<vector<handle_t>> consolidate_bridges(const HandleGraph& graph,
                                             const vector<handle_t>& bridges) {
    
    // return a bool if the node has degree 1 in that direction, and if so also
    // the neighbor node
    auto degree_one_neighbor = [&](handle_t n, bool left_side) {
        int count = 0;
        handle_t nbr;
        graph.follow_edges(n, left_side, [&](const handle_t& next) {
            ++count;
            nbr = next;
            return count <= 1;
        });
        return make_pair(count == 1, nbr);
    };
        
    unordered_set<handle_t> remaining(bridges.begin(), bridges.end());
    
    vector<vector<handle_t>> return_val;
    
    for (handle_t node : bridges) {
        
        if (!remaining.count(node)) {
            continue;
        }
        
        // we will extend this node into a consolidated bridge
        
        return_val.emplace_back(1, node);
        auto& consolidated_bridge = return_val.back();
        remaining.erase(node);
        
        // try to extend the bridge in both directions
        for (bool to_left : {true, false}) {
            
            while (true) {
                
                bool degree_one;
                handle_t nbr;
                tie(degree_one, nbr) = degree_one_neighbor(consolidated_bridge.back(), to_left);
                if (!degree_one || !remaining.count(graph.forward(nbr))) {
                    // there is not a single neighbor, or else it is not a bridge
                    break;
                }
                if (!degree_one_neighbor(nbr, !to_left).first) {
                    // the neighbor does not have degree 1 on this side
                    break;
                }
                consolidated_bridge.push_back(nbr);
                remaining.erase(nbr);
            }
            
            // when we build leftwards, we need to reverse it to keep the
            // order consistent
            if (to_left) {
                reverse(consolidated_bridge.begin(), consolidated_bridge.end());
            }
        }
    }
    
    assert(remaining.empty());
    return return_val;
}

// FIXME: this will not find any bridge-free connected components...

void for_each_bridge_component(const HandleGraph& graph,
                               const vector<vector<handle_t>>& bridges,
                               const function<void(const HandleGraph&,
                                                   const vector<pair<size_t, bool>>&)>& f) {
    
    
    unordered_map<handle_t, pair<size_t, bool>> bridge_ends;
    unordered_set<handle_t> visited;
    for (size_t i = 0; i < bridges.size(); ++i) {
        const auto& bridge = bridges[i];
        bridge_ends[bridge.back()] = make_pair(i, false);
        bridge_ends[graph.flip(bridge.front())] = make_pair(i, true);
        
        // we pre-mark the sides we don't want to traverse as visited
        for (size_t j = 0; j < bridge.size(); ++j) {
            if (j != 0) {
                visited.insert(graph.flip(bridge[j]));
            }
            if (j + 1 != bridge.size()) {
                visited.insert(bridge[j]);
            }
        }
    }
    
    
    graph.for_each_handle([&](const handle_t& node) {
        for (bool left_side : {false, true}) {
            
            handle_t seed = left_side ? graph.flip(node) : node;
            if (visited.count(seed)) {
                // we've already traversed this side's component
                continue;
            }
            
            visited.insert(seed);
            
            vector<pair<size_t, bool>> adjacent_bridges;
            unordered_set<nid_t> node_ids{graph.get_id(seed)};
            vector<handle_t> stack(1, seed);
            
            while (!stack.empty()) {
                
                handle_t node = stack.back();
                stack.pop_back();
                
                auto it = bridge_ends.find(node);
                if (it != bridge_ends.end()) {
                    // we hit a bridge, don't cross it, but remember which
                    // bridge it was
                    adjacent_bridges.push_back(it->second);
                }
                else {
                    // this is not a bridge, we can look at the other side
                    handle_t across_node = graph.flip(node);
                    if (!visited.count(across_node)) {
                        stack.push_back(across_node);
                        visited.insert(across_node);
                    }
                }
                
                // we can always look across edges
                graph.follow_edges(node, false, [&](const handle_t& next) {
                    handle_t adjacent = graph.flip(next);
                    if (!visited.count(adjacent)) {
                        node_ids.insert(graph.get_id(adjacent));
                        visited.insert(adjacent);
                        stack.push_back(adjacent);
                    }
                });
            }
            
            // make a component subgraph and execute the lambda
            SubgraphOverlay bridge_component(&graph, &node_ids);
            f(bridge_component, adjacent_bridges);
        }
    });
}

}
