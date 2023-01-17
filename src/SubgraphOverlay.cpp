#include <atomic>
#include "SubgraphOverlay.hpp"

namespace gfase {

SubgraphOverlay::SubgraphOverlay(const HandleGraph& backing, const unordered_set<nid_t>& node_subset) :
    backing_graph(&backing),
    node_subset(&node_subset) {
    if (!node_subset->empty()) {
        auto minmax_nodes = std::minmax_element(node_subset->begin(), node_subset->end());
        min_node = *minmax_nodes.first;
        max_node = *minmax_nodes.second;
    } 
}

SubgraphOverlay::~SubgraphOverlay() {
    
}

bool SubgraphOverlay::has_node(nid_t node_id) const {
    return node_subset->count(node_id);
}
   
handle_t SubgraphOverlay::get_handle(const nid_t& node_id, bool is_reverse) const {
    if (has_node(node_id)) {
        return backing_graph->get_handle(node_id, is_reverse);
    } else {
        throw runtime_error("Node " + std::to_string(node_id) + " not in subgraph overlay");
    }
}
    
nid_t SubgraphOverlay::get_id(const handle_t& handle) const {
    return backing_graph->get_id(handle);
}
    
bool SubgraphOverlay::get_is_reverse(const handle_t& handle) const {
    return backing_graph->get_is_reverse(handle);
}

handle_t SubgraphOverlay::flip(const handle_t& handle) const {
    return backing_graph->flip(handle);
}
    
size_t SubgraphOverlay::get_length(const handle_t& handle) const {
    return backing_graph->get_length(handle);
}

std::string SubgraphOverlay::get_sequence(const handle_t& handle) const {
    return backing_graph->get_sequence(handle);
}
    
size_t SubgraphOverlay::get_node_count() const {
    return node_subset->size();
}

nid_t SubgraphOverlay::min_node_id() const {
    return min_node;
}
    
nid_t SubgraphOverlay::max_node_id() const {
    return max_node;
}

bool SubgraphOverlay::follow_edges_impl(const handle_t& handle, bool go_left, const std::function<bool(const handle_t&)>& iteratee) const {
    
    std::function<bool(const handle_t&)> subgraph_iteratee = [&](const handle_t& handle) {
        if (has_node(backing_graph->get_id(handle))) {
            if (!iteratee(handle)) {
                return false;
            }
        }
        return true;
    };
    if (has_node(backing_graph->get_id(handle))) {
        return backing_graph->follow_edges(handle, go_left, subgraph_iteratee);
    }
    return true;
}
    
bool SubgraphOverlay::for_each_handle_impl(const std::function<bool(const handle_t&)>& iteratee, bool parallel) const {

    if (!parallel) {
        bool keep_going = true;
        for (auto node_it = node_subset->begin(); keep_going && node_it != node_subset->end(); ++node_it) {
            keep_going = iteratee(get_handle(*node_it, false));
        }
        return keep_going;
    } else {
        // copy them into something easy to iterate with omp
        vector<nid_t> node_vec(node_subset->begin(), node_subset->end());
        std::atomic<bool> keep_going(true);
#pragma omp parallel for
        for (size_t i = 0; i < node_vec.size(); ++i) {
            keep_going = keep_going && iteratee(backing_graph->get_handle(node_vec[i]));
        }
        return keep_going;
    }
}

}
