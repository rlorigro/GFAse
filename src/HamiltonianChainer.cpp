#include "HamiltonianChainer.hpp"

#include "HamiltonianPath.hpp"
#include "Bridges.hpp"
#include "graph_utility.hpp"

#include "handlegraph/util.hpp"
#include "handlegraph/types.hpp"

using handlegraph::as_handle;
using handlegraph::handle_t;
using handlegraph::path_handle_t;

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <cassert>
#include <deque>
#include <array>
#include <algorithm>

using std::vector;
using std::unordered_set;
using std::unordered_map;
using std::set;
using std::move;
using std::deque;
using std::array;
using std::reverse;

namespace gfase {

/*
 * A 0-, 1-, or 2-sided component with 2 fully or partially walked alleles
 */
struct ChainableComponent {
    ChainableComponent() = default;
    ~ChainableComponent() = default;
    
    // we assume that the left side is assigned first (i.e. there can be a left
    // side and no right side, but not the reverse)
    // the sides point are oriented to point into the component
    bool has_left_side = false;
    bool has_right_side = false;
    handle_t left_side = as_handle(-1);
    handle_t right_side = as_handle(-1);
    
    // the longest certain alleles we can walk inward from the left boundaries.
    // if there are no boundaries, these start from an arbitrary node.
    // includes the boundary node(s) when they are present
    array<vector<handle_t>, 2> allele_from_left;
    
    // the longest certain alleles we can walk inward from the right boundaries.
    // left blank if there is no right boundary, or if the left-side allele spans
    // the entire component
    array<vector<handle_t>, 2> allele_from_right;
    
};


void HamiltonianChainer::generate_chain_paths(MutablePathDeletableHandleGraph& graph,
                                              const MultiContactGraph& contact_graph) const {
    
    
    /*
     * Part 1: find bridges
     */
    
    // find bridges that are not assigned to one or the other haplotype
    vector<handle_t> nonhaploid_bridges;
    for_each_connected_component_subgraph(graph, [&](const HandleGraph& component) {
        for (handle_t bridge : bridge_nodes(component)) {
            if (!contact_graph.has_node(graph.get_id(bridge)) ||
                !contact_graph.has_alt(graph.get_id(bridge))) {
                nonhaploid_bridges.push_back(bridge);
            }
        }
    });
    
    // merge into unipath bridges
    auto unipath_bridges = consolidate_bridges(graph, nonhaploid_bridges);
    
    // and we can use this to look up bridges by their component-facing boundaries
    unordered_map<handle_t, size_t> boundary_to_bridge;
    boundary_to_bridge.reserve(2 * unipath_bridges.size());
    for (size_t i = 0; i < unipath_bridges.size(); ++i) {
        const auto& unipath_bridge = unipath_bridges[i];
        boundary_to_bridge[graph.flip(unipath_bridge.front())] = i;
        boundary_to_bridge[unipath_bridge.back()] = i;
    }
    
    /*
     * Part 2: identify walks through bridge components with bridge degree <= 2 and
     * non-bridge simple bubbles
     */
    
    // we'll try to fill these out for each bridge component
    unordered_map<handle_t, size_t> boundary_to_chain_link;
    vector<ChainableComponent> chain_links;
    
    vector<vector<handle_t>> bubble_unipath_boundaries;
    
    for_each_bridge_component(graph, unipath_bridges,
                              [&](const HandleGraph& bridge_component,
                                  const vector<pair<size_t, bool>>& incident_bridges) {
        
        
        bool found_allele_success = false;
        if (incident_bridges.size() <= 2) {
            
            // this bridge component is a potentially phaseable unit in itself
            
            // check that the bridge component looks like it's capturing two allelic sequences
            unordered_set<nid_t> phase_0_nodes, phase_1_nodes;
            bool all_phasable = bridge_component.for_each_handle([&](const handle_t& handle) {
                
                nid_t node_id = bridge_component.get_id(handle);
                if (phase_0_nodes.count(node_id) || phase_1_nodes.count(node_id)) {
                    // we already processed this node from another node's alt component
                    return true;
                }
                if (contact_graph.has_node(node_id) && contact_graph.has_alt(node_id)) {
                    alt_component_t alt_component;
                    contact_graph.get_alt_component(node_id, false, alt_component);
                    for (auto alt_set : {&alt_component.first, &alt_component.second}) {
                        for (auto member_id : *alt_set) {
                            if (!bridge_component.has_node(member_id)) {
                                // a member of this alt set is not found in the bridge component.
                                // this makes it less likely that the bridge component represents
                                // a pair of allelic sequences, but we'll still allow it as long as
                                // the missing sequences are isolated nodes
                                // TODO: this condition is motivated by patterns we've seen in
                                // verkko graphs, but it could stand to be a bit more principled
                                handle_t missing_node = graph.get_handle(member_id);
                                bool no_edges = graph.follow_edges(missing_node, false, [&](const handle_t& null) {
                                    return false;
                                });
                                no_edges = no_edges && graph.follow_edges(missing_node, true, [&](const handle_t& null) {
                                    return false;
                                });
                                if (!no_edges) {
                                    // this bridge component is not phasable
                                    return false;
                                }
                            }
                        }
                    }
                    // record the phase of the two alt sets (which also marks them as having been processed)
                    bool order_swapped = (contact_graph.get_partition(*alt_component.first.begin()) > 0);
                    for (auto alt_allele_id : (order_swapped ? alt_component.second : alt_component.first)) {
                        if (bridge_component.has_node(alt_allele_id)) {
                            phase_0_nodes.insert(alt_allele_id);
                        }
                    }
                    for (auto alt_allele_id : (order_swapped ? alt_component.first : alt_component.second)) {
                        if (bridge_component.has_node(alt_allele_id)) {
                            phase_1_nodes.insert(alt_allele_id);
                        }
                    }
                }
                return true;
            });
            
            if (!all_phasable) {
                // the nodes had alts that are outside this bridge component
                return;
            }
            
            unordered_set<handle_t> start, end;
            if (incident_bridges.size() >= 1) {
                // we'll start at the inward side of the bridge, facing into to the component
                const auto& start_bridge = unipath_bridges[incident_bridges[0].first];
                if (incident_bridges[0].second) {
                    start.insert(graph.flip(start_bridge.front()));
                }
                else {
                    start.insert(start_bridge.back());
                }
            }
            if (incident_bridges.size() == 2) {
                // we'll end at the inward side of the bridge, but facing outward
                const auto& end_bridge = unipath_bridges[incident_bridges[1].first];
                if (incident_bridges[1].second) {
                    end.insert(end_bridge.front());
                }
                else {
                    end.insert(graph.flip(end_bridge.back()));
                }
            }
            
            bool resolved_hamiltonian_0;
            auto phase_0_walks = generate_allelic_semiwalks(bridge_component,
                                                            phase_0_nodes,
                                                            phase_1_nodes,
                                                            start, end,
                                                            resolved_hamiltonian_0);
            
            bool resolved_hamiltonian_1;
            auto phase_1_walks = generate_allelic_semiwalks(bridge_component,
                                                            phase_1_nodes,
                                                            phase_0_nodes,
                                                            start, end,
                                                            resolved_hamiltonian_1);
            
            // we'll consider this bridge component fully solved (even without unique alleles)
            // if we found a hamiltonian path for either of the phases
            found_allele_success = (resolved_hamiltonian_0 || resolved_hamiltonian_1);
            
            // record the result of the allele identification
            chain_links.emplace_back();
            auto& link = chain_links.back();
            if (!start.empty()) {
                link.has_left_side = true;
                link.left_side = *start.begin();
                boundary_to_chain_link[link.left_side] = chain_links.size() - 1;
            }
            if (!end.empty()) {
                link.has_right_side = true;
                link.right_side = bridge_component.flip(*end.begin());
                boundary_to_chain_link[link.right_side] = chain_links.size() - 1;
            }
            link.allele_from_left[0] = move(phase_0_walks.first);
            link.allele_from_right[0] = move(phase_0_walks.second);
            link.allele_from_left[1] = move(phase_1_walks.first);
            link.allele_from_right[1] = move(phase_1_walks.second);
        }
        
        if (incident_bridges.size() > 2 || !found_allele_success) {
            // this bridge component has bridge degree > 2 or else the hamiltonian algorithm failed on both
            // alleles. we could maybe find phaseable ubbles inside it using a rigid topological motif criterion
            
            unordered_set<handle_t> processed_sides;
            bridge_component.for_each_handle([&](const handle_t& handle) {
                nid_t node_id = bridge_component.get_id(handle);
                if (contact_graph.has_node(node_id) && contact_graph.has_alt(node_id)) {
                    // we don't want to find haploid allelic sequences, we want boundaries of bubbles
                    return;
                }
                for (auto side : {handle, bridge_component.flip(handle)}) {
                    if (processed_sides.count(side)) {
                        continue;
                    }
                    // check if this is the boundary of a bubble
                    processed_sides.insert(side);
                    
                    vector<handle_t> neighbors;
                    bridge_component.follow_edges(side, false, [&](const handle_t& nbr) {
                        neighbors.push_back(nbr);
                    });
                    if (neighbors.size() != 2) {
                        // these can't be the two sides of a bubble
                        continue;
                    }
                    array<int, 2> neighbor_rev_count{0, 0};
                    for (int hap : {0, 1}) {
                        bridge_component.follow_edges(neighbors[hap], true, [&](const handle_t& prev) {
                            ++neighbor_rev_count[hap];
                            return neighbor_rev_count[hap] < 2;
                        });
                    }
                    if (neighbor_rev_count[0] != 1 || neighbor_rev_count[1] != 1) {
                        // you can reach the two nodes from other nodes than the boundary
                        continue;
                    }
                    array<handle_t, 2> next_neighbor{as_handle(-1), as_handle(-1)};
                    bool deg_less_than_2 = true;
                    for (int hap : {0, 1}) {
                        deg_less_than_2 = deg_less_than_2 && bridge_component.follow_edges(neighbors[hap], false,
                                                                                           [&](const handle_t& next) {
                            if (next_neighbor[hap] != as_handle(-1)) {
                                return false;
                            }
                            else {
                                next_neighbor[hap] = next;
                                return true;
                            }
                        });
                    }
                    if (next_neighbor[0] == as_handle(-1) || next_neighbor[0] != next_neighbor[1] || !deg_less_than_2) {
                        // they don't meet together at the following node
                        continue;
                    }
                    size_t next_neighbor_rev_count = 0;
                    bridge_component.follow_edges(next_neighbor[0], true, [&](const handle_t& prev) {
                        ++next_neighbor_rev_count;
                        return next_neighbor_rev_count < 3;
                    });
                    if (next_neighbor_rev_count > 2) {
                        // other nodes also meet together at the following node
                        continue;
                    }
                    
                    if (!contact_graph.has_node(bridge_component.get_id(neighbors[0])) ||
                        !contact_graph.has_node(bridge_component.get_id(neighbors[1]))) {
                        // these can't have an assigned phase
                        continue;
                    }
                    alt_component_t alt_component;
                    contact_graph.get_alt_component(bridge_component.get_id(neighbors[0]), false, alt_component);
                    if (alt_component.first.size() != 1 || alt_component.second.size() != 1) {
                        // these have other homology partners outside of this motif
                        continue;
                    }
                    if ((*alt_component.first.begin() != bridge_component.get_id(neighbors[0]) ||
                         *alt_component.second.begin() != bridge_component.get_id(neighbors[1])) &&
                        (*alt_component.second.begin() != bridge_component.get_id(neighbors[0]) ||
                         *alt_component.first.begin() != bridge_component.get_id(neighbors[1]))) {
                        // they aren't homology partners with each other
                        continue;
                    }
                    if (contact_graph.has_node(bridge_component.get_id(next_neighbor[0])) &&
                        contact_graph.has_alt(bridge_component.get_id(next_neighbor[0]))) {
                        // the other boundary is phased (i.e. isn't diploid)
                        continue;
                    }
                    
                    // we've fully checked the local topology, this looks like a phased bubble
                    
                    chain_links.emplace_back();
                    auto& link = chain_links.back();
                    link.has_left_side = true;
                    link.left_side = side;
                    link.has_right_side = true;
                    link.right_side = bridge_component.flip(next_neighbor[0]);
                    
                    // form the alleles
                    vector<handle_t> allele_0{side, neighbors[0], next_neighbor[0]};
                    vector<handle_t> allele_1{side, neighbors[1], next_neighbor[0]};
                    bool order_swapped = (contact_graph.get_partition(bridge_component.get_id(neighbors[0])) > 0);
                    link.allele_from_left[0] = (order_swapped ? move(allele_1) : move(allele_0));
                    link.allele_from_left[1] = (order_swapped ? move(allele_0) : move(allele_1));
                    
                    // walk out the full unipath boundary (i.e. not just the inward-facing node)
                    bubble_unipath_boundaries.push_back(walk_diploid_unipath(bridge_component, contact_graph,
                                                                             link.left_side, true));
                    bubble_unipath_boundaries.push_back(walk_diploid_unipath(bridge_component, contact_graph,
                                                                             link.right_side, false));
                    
                    processed_sides.insert(side);
                    for (int i = 0; i < 2; ++i) {
                        processed_sides.insert(neighbors[i]);
                        processed_sides.insert(bridge_component.flip(neighbors[i]));
                    }
                    processed_sides.insert(bridge_component.flip(next_neighbor[0]));
                }
            });
        }
    });
    
    // move the unipath boundaries from the bubbles to the same list as the bridges (even though
    // they're technically not bridges)
    // note: we don't need to worry about re-identifuying a bridge because if a bridge satisfied
    // the rigid topological bubble criterion, then it would have been solved by the hamiltonian
    for (auto& unipath_bridge : bubble_unipath_boundaries) {
        boundary_to_bridge[graph.flip(unipath_bridge.front())] = unipath_bridges.size();
        boundary_to_bridge[unipath_bridge.back()] = unipath_bridges.size();
        unipath_bridges.emplace_back(move(unipath_bridge));
    }
    
    /*
     * Part 3: find bridges
     */
    
    array<int, 2> next_path_id{0, 0};
    auto get_next_path_name = [&](int hap) {
        return string("gfase_hap_" + to_string(next_path_id[hap]++) +"_" + to_string(hap));
    };
    
    vector<bool> is_chained(chain_links.size(), false);
    
    // add the next bridge and chain link to the growing chain
    // TODO: if i split this into the bridge and the component walk, i could avoid doing the same
    // bridge traversing logic on both haplotypes...
    auto add_next_component = [&](int haplotype,
                                  bool extending_backward,
                                  bool& left_to_right, // which direction is "forward" relative to the path
                                  path_handle_t& curr_path,
                                  int64_t& curr_link_idx,
                                  bool& keep_going,
                                  bool& found_cycle) {
        
        auto& curr_link = chain_links[curr_link_idx];
        bool go_to_left = (extending_backward == left_to_right);
        if ((!curr_link.has_left_side && go_to_left) || (!curr_link.has_right_side && !go_to_left)) {
            // TODO: should i actually mark that we don't keep going after adding the allele?
            keep_going = false;
            return;
        }
        
        // add the nodes of the bridge
        handle_t adjacent_side = go_to_left ? curr_link.left_side : curr_link.right_side;
        size_t bridge_idx = boundary_to_bridge.at(adjacent_side);
        auto& unipath_bridge = unipath_bridges[bridge_idx];
        handle_t final_outward; // the last node of the bridge, facing into the next component
        if (unipath_bridge.back() == adjacent_side) {
            // we're traversing the bridge backwards
            for (int64_t j = unipath_bridges.size() - 2; j >= 0; --j) {
                if (extending_backward) {
                    graph.prepend_step(curr_path, unipath_bridge[j]);
                }
                else {
                    graph.append_step(curr_path, graph.flip(unipath_bridge[j]));
                }
            }
            final_outward = graph.flip(unipath_bridge.front());
        }
        else {
            // we're traversing the bridge forwards
            assert(adjacent_side == graph.flip(unipath_bridge.front()));
            for (size_t j = 1; j < unipath_bridges.size(); ++j) {
                if (extending_backward) {
                    graph.prepend_step(curr_path, graph.flip(unipath_bridge[j]));
                }
                else {
                    graph.append_step(curr_path, unipath_bridge[j]);
                }
            }
            final_outward = unipath_bridge.back();
        }
        
        // gather information about the next chain link and its orientation
        curr_link_idx = boundary_to_chain_link.at(final_outward);
        if (is_chained[curr_link_idx]) {
            // we've previously walked this chain (probably it's in a cycle) and don't want
            // to do it again
            // FIXME: should we make the path be cyclic then?
            keep_going = false;
            found_cycle = true;
            return;
        }
        auto& next_link = chain_links[curr_link_idx];
        assert(final_outward == next_link.left_side || final_outward == next_link.right_side);
        bool enter_left = final_outward == next_link.left_side;
        left_to_right = enter_left != extending_backward;
        
        // handles the logic of adding either the left or right allele onto the path, optionally
        // checks consistency of whether
        auto add_allele = [&](const vector<handle_t>& allele, handle_t* check_inward) {
            if (enter_left) {
                // we traverse this component left to right
                if (check_inward) {
                    assert(*check_inward == (left_to_right ? allele.front() : graph.flip(allele.front())));
                }
                for (size_t i = 1; i < allele.size(); ++i) {
                    if (extending_backward) {
                        graph.prepend_step(curr_path, graph.flip(allele[i]));
                    }
                    else {
                        graph.append_step(curr_path, allele[i]);
                    }
                }
            }
            else {
                // we traverse this component right to left
                if (check_inward) {
                    assert(*check_inward == (left_to_right ? allele.back() : graph.flip(allele.back())));
                }
                for (int64_t i = allele.size() - 2; i >= 0; --i) {
                    if (extending_backward) {
                        graph.prepend_step(curr_path, allele[i]);
                    }
                    else {
                        graph.append_step(curr_path, graph.flip(allele[i]));
                    }
                }
            }
        };
        
        // add the alleles to the path
        if (next_link.allele_from_right[haplotype].empty()) {
            // the left walk goes across the entire link, or else there is no right side
            add_allele(next_link.allele_from_left[haplotype], &final_outward);
            keep_going = next_link.has_right_side;
        }
        else {
            // the walk is broken in the middle
            auto& left_allele = next_link.allele_from_left[haplotype];
            auto& right_allele = next_link.allele_from_right[haplotype];
            if (enter_left) {
                add_allele(left_allele, &final_outward);
                // end the haplotype path and start a new one
                curr_path = graph.create_path_handle(get_next_path_name(haplotype));
                add_allele(right_allele, nullptr);
            }
            else {
                add_allele(right_allele, &final_outward);
                // end the haplotype path and start a new one
                curr_path = graph.create_path_handle(get_next_path_name(haplotype));
                add_allele(left_allele, nullptr);
            }
        }
        
        is_chained[curr_link_idx] = true;
    };
    
    for (size_t i = 0; i < chain_links.size(); ++i) {
        if (is_chained[i]) {
            continue;
        }
        
        // make paths to extend into haplotypes
        array<path_handle_t, 2> init_paths;
        for (int hap : {0, 1}) {
            init_paths[hap] = graph.create_path_handle(get_next_path_name(hap));
        }
        array<path_handle_t, 2> curr_paths = init_paths;
        auto& init_link = chain_links[i];
        for (int hap : {0, 1}) {
            for (handle_t step : init_link.allele_from_left[hap]) {
                graph.append_step(curr_paths[hap], step);
            }
        }
        // TODO: it's a bit silly that i maintain these variables separately for both haplotypes,
        // but i don't want it to get set during one haplotype iteration and mess up the next
        array<bool, 2> link_is_left_to_right{true, true};
        array<int64_t, 2> link_index{(int64_t) i, (int64_t) i};
        bool keep_going = init_link.has_left_side;
        bool found_cycle = false;
        
        // start by extending to the left
        while (keep_going) {
            for (int hap : {0, 1}) {
                add_next_component(hap, true, link_is_left_to_right[hap],
                                   curr_paths[hap], link_index[hap], keep_going, found_cycle);
            }
        }
        
        if (found_cycle) {
            // we looped back around to the same link, but we might not have included the
            // right side of this link
            for (int hap : {0, 1}) {
                auto& allele = init_link.allele_from_right[hap];
                if (!allele.empty()) {
                    for (int64_t j = allele.size() - 2; j >= 0; --j) {
                        graph.prepend_step(curr_paths[hap], allele[j]);
                    }
                }
            }
        }
        else if (init_link.has_right_side) {
            // now we extend to the right
            
            // initialize the rightward traversal
            keep_going = true;
            for (int hap : {0, 1}) {
                link_is_left_to_right[hap] = true;
                link_index[hap] = i;
                if (init_link.allele_from_right[hap].empty()) {
                    // we can continue the initial allele
                    curr_paths[hap] = init_paths[hap];
                }
                else {
                    // we have to start a new path
                    curr_paths[hap] = graph.create_path_handle(get_next_path_name(hap));
                }
                for (handle_t step : init_link.allele_from_right[hap]) {
                    graph.append_step(curr_paths[hap], step);
                }
            }
            
            // add links to the right
            while (keep_going) {
                for (int hap : {0, 1}) {
                    add_next_component(hap, false, link_is_left_to_right[hap],
                                       curr_paths[hap], link_index[hap], keep_going, found_cycle);
                }
            }
        }
        is_chained[i] = true;
    }
}

pair<vector<handle_t>, vector<handle_t>>
HamiltonianChainer::generate_allelic_semiwalks(const HandleGraph& graph,
                                               const unordered_set<nid_t>& in_phase_nodes,
                                               const unordered_set<nid_t>& out_phase_nodes,
                                               const unordered_set<handle_t>& starts,
                                               const unordered_set<handle_t>& ends,
                                               bool& resolved_hamiltonian) const {
    
    pair<vector<handle_t>, vector<handle_t>> return_val;
    
    // for reverse iteration (if we need to do it)
    unordered_set<handle_t> rev_starts, rev_ends;
    for (auto end : ends) {
        rev_starts.insert(graph.flip(end));
    }
    for (auto start : starts) {
        rev_ends.insert(graph.flip(start));
    }
    
    auto hamiltonian = find_hamiltonian_path(graph,
                                             in_phase_nodes,
                                             out_phase_nodes,
                                             starts, ends, hamiltonian_max_iters);
    
    if (hamiltonian.is_solved && !hamiltonian.hamiltonian_path.empty()) {
        // this phase can be walked out as a hamiltonian path
        resolved_hamiltonian = true;
        
        // is the hamiltonian unique?
        bool fully_unique = (hamiltonian.unique_prefix.size() == hamiltonian.hamiltonian_path.size());
        return_val.first = move(hamiltonian.unique_prefix);
        
        if (!fully_unique && !ends.empty()) {
            // we'll also try to get the unique parts of the other end
            auto rev_hamiltonian = find_hamiltonian_path(graph,
                                                         in_phase_nodes,
                                                         out_phase_nodes,
                                                         rev_starts, rev_ends, hamiltonian_max_iters);
            
            // TODO: does solvability within the max iters guarantee it for the other side?
            // should I additionally fall back to the unambiguous walk if this fails?
            if (rev_hamiltonian.is_solved) {
                // reverse the prefix (to turn it into a suffix) and add it to the return value
                auto& prefix = rev_hamiltonian.unique_prefix;
                for (auto it = prefix.rbegin(); it != prefix.rend(); ++it) {
                    return_val.second.push_back(graph.flip(*it));
                }
            }
        }
        
        // TODO: should I treat it as ambiguous if there are cycles that involve the chosen path
        // (in which case the correct walk might not be Hamiltonian)
    }
    else {
        // fall back on completely unambiguous walks if we can't find a full walk
        resolved_hamiltonian = false;
        if (starts.size() == 1) {
            return_val.first = max_unambiguous_path(graph, *starts.begin(), in_phase_nodes);
        }
        if (ends.size() == 1) {
            auto rev_suffix = max_unambiguous_path(graph, *rev_starts.begin(), in_phase_nodes);
            for (auto it = rev_suffix.rbegin(); it != rev_suffix.rend(); ++it) {
                return_val.second.push_back(graph.flip(*it));
            }
        }
    }
    return return_val;
}


vector<handle_t> HamiltonianChainer::max_unambiguous_path(const HandleGraph& graph, handle_t start,
                                                          const unordered_set<nid_t>& allowed_nodes) const {
    
    vector<handle_t> walk(1, start);
    bool unambiguous = true;
    while (unambiguous) {
        handle_t to_add = as_handle(-1);
        size_t num_neighbors = 0;
        unambiguous = graph.follow_edges(walk.back(), false, [&](const handle_t& next) {
            ++num_neighbors;
            to_add = next;
            return num_neighbors == 1 && allowed_nodes.count(graph.get_id(next));
        });
        if (num_neighbors != 0 && unambiguous) {
            walk.push_back(to_add);
        }
    }
    return walk;
}

vector<handle_t> HamiltonianChainer::walk_diploid_unipath(const HandleGraph& graph, const MultiContactGraph& contact_graph,
                                                          handle_t handle, bool go_left) const {
    vector<handle_t> unipath(1, handle);
    while (true) {
        int degree = 0;
        handle_t neighbor = as_handle(-1);
        graph.follow_edges(unipath.back(), go_left, [&](const handle_t& next) {
            ++degree;
            neighbor = next;
            return degree < 2;
        });
        if (degree != 1) {
            break;
        }
        degree = 0;
        graph.follow_edges(neighbor, !go_left, [&](const handle_t& prev) {
            ++degree;
            return degree < 2;
        });
        if (degree != 1) {
            break;
        }
        if (contact_graph.has_node(graph.get_id(neighbor)) &&
            contact_graph.has_alt(graph.get_id(neighbor))) {
            break;
        }
        unipath.push_back(neighbor);
    }
    if (go_left) {
        reverse(unipath.begin(), unipath.end());
    }
    return unipath;
}

}
