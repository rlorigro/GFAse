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
#include <iostream>
#include <stdexcept>

using std::vector;
using std::unordered_set;
using std::unordered_map;
using std::set;
using std::move;
using std::deque;
using std::array;
using std::reverse;
using std::cerr;
using std::endl;
using std::runtime_error;

namespace gfase {

static const bool debug = false;

/*
 * A 0-, 1-, or 2-sided component with 2 fully or partially walked alleles
 */
struct ChainableComponent {
    ChainableComponent() = default;
    ~ChainableComponent() = default;
    
    // we assume that the left side is assigned first (i.e. there can be a left
    // side and no right side, but not the reverse)
    // the sides are oriented to point into the component
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

bool HamiltonianChainer::has_phase_chain(const string& name) const {
    if (!path_graph || !path_graph->has_path(name)) {
        return false;
    }
    path_handle_t path_handle = path_graph->get_path_handle(name);
    return (phase_paths[0].count(path_handle) || phase_paths[1].count(path_handle));
}

int8_t HamiltonianChainer::get_partition(const string& name) const {
    if (!path_graph || !path_graph->has_path(name)) {
        return 0;
    }
    path_handle_t path_handle = path_graph->get_path_handle(name);
    if (phase_paths[0].count(path_handle)) {
        return -1;
    }
    else if (phase_paths[1].count(path_handle)) {
        return 1;
    }
    else {
        return 0;
    }
}

void HamiltonianChainer::generate_chain_paths(MutablePathDeletableHandleGraph& graph,
                                              const IncrementalIdMap<string>& id_map,
                                              const MultiContactGraph& contact_graph) {
    
    
    // TODO: this is ugly, but it saves me from having to re-creating indexes that are mostly identical
    // to the path handle interface
    if (path_graph) {
        throw runtime_error("ERROR: attempted to generate chain paths twice with the same HamiltonianChainer");
    }
    path_graph = &graph;
    
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
    if (debug) {
        cerr << "identified unipath bridges:" << endl;
        for (const auto& bridge : unipath_bridges) {
            cerr << "\t";
            for (size_t i = 0; i < bridge.size(); ++i) {
                if (i != 0) {
                    cerr << ", ";
                }
                cerr << graph.get_id(bridge[i]) << (graph.get_is_reverse(bridge[i]) ? "-" : "+") << "(" << id_map.get_name(graph.get_id(bridge[i])) << ")";
            }
            cerr << endl;
        }
    }
    
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
     * non-bridge-bordered simple bubbles
     */
    
    // we'll try to fill these out for each bridge component
    unordered_map<handle_t, size_t> boundary_to_chain_link;
    vector<ChainableComponent> chain_links;
    
    vector<vector<handle_t>> bubble_unipath_boundaries;
    
    for_each_bridge_component(graph, unipath_bridges,
                              [&](const HandleGraph& bridge_component,
                                  const vector<pair<size_t, bool>>& incident_bridges) {
        
        if (debug) {
            cerr << "phasing in component bordering bridges:" << endl;
            for (const auto& bridge_idx : incident_bridges) {
                cerr << "\trev? " << bridge_idx.second << ": ";
                const auto& bridge = unipath_bridges[bridge_idx.first];
                for (size_t i = 0; i < bridge.size(); ++i) {
                    if (i != 0) {
                        cerr << ", ";
                    }
                    cerr << graph.get_id(bridge[i]) << (graph.get_is_reverse(bridge[i]) ? "-" : "+");
                }
                cerr << endl;
            }
        }
        
        bool found_allele_success = false;
        if (incident_bridges.size() <= 2) {
            if (debug) {
                cerr << "bridge degree is " << incident_bridges.size() << ", attempting to find hamiltonian alleles" << endl;
            }
            
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
                                                            id_map,
                                                            phase_0_nodes,
                                                            phase_1_nodes,
                                                            start, end,
                                                            resolved_hamiltonian_0);
            
            bool resolved_hamiltonian_1;
            auto phase_1_walks = generate_allelic_semiwalks(bridge_component,
                                                            id_map,
                                                            phase_1_nodes,
                                                            phase_0_nodes,
                                                            start, end,
                                                            resolved_hamiltonian_1);
            
            // we'll consider this bridge component fully solved (even without unique alleles)
            // if we found a hamiltonian path for either of the phases
            found_allele_success = (resolved_hamiltonian_0 || resolved_hamiltonian_1);
            
            if (found_allele_success) {
                // TODO: is this the right logic? should i ever add the partial alleles even if
                // the hamiltonian isn't successful? it's hard to know when to fall back on the
                // smaller simple bubbles contained in the component...
                
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
                
                if (debug) {
                    cerr << "succeeded in finding hamiltonian allele(s), added a chainable component" << endl;
                    for (auto left : {true, false}) {
                        cerr << "from " << (left ? "left" : "right") << ":" << endl;
                        auto alleles = left ? link.allele_from_left : link.allele_from_right;
                        for (auto allele : alleles) {
                            for (auto handle : allele) {
                                cerr << " " << bridge_component.get_id(handle) << "(" << id_map.get_name(bridge_component.get_id(handle)) << ")" << (bridge_component.get_is_reverse(handle) ? "-" : "+");
                            }
                            cerr << endl;
                        }
                    }
                }
            }
        }
        
        if (incident_bridges.size() > 2 || !found_allele_success) {
            // this bridge component has bridge degree > 2 or else the hamiltonian algorithm failed on both
            // alleles. we could maybe find phaseable bubbles inside it using a rigid topological motif criterion
            
            if (debug) {
                cerr << "no hamiltonian alleles are possible, attempting to find chainable simple bubbles" << endl;
            }
            
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
                    
                    if (debug) {
                        cerr << "found a chainable simple bubble consisting of alleles:" << endl;
                        for (auto allele : link.allele_from_left) {
                            for (auto handle : allele) {
                                cerr << " " << bridge_component.get_id(handle) << (bridge_component.get_is_reverse(handle) ? "-" : "+");
                            }
                            cerr << endl;
                        }
                    }
                    
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
    
    if (debug) {
        cerr << "found " << chain_links.size() << " chainable components" << endl;
    }
    
    /*
     * Part 3: form chains
     */
    
    auto get_next_path = [&](int hap) {
        path_handle_t path_handle = graph.create_path_handle(phase_path_name(hap, next_path_ids[hap]++));
        if (debug) {
            cerr << "generating new path named " << graph.get_path_name(path_handle) << endl;
        }
        phase_paths[hap].insert(path_handle);
        return path_handle;
    };
    
    vector<bool> is_chained(chain_links.size(), false);
    
    // function to add the next bridge to both chains and advance the tracking variables to the next link
    auto add_next_bridge = [&](bool extending_backward,
                               bool& left_to_right,
                               const array<path_handle_t, 2>& curr_paths,
                               int64_t& curr_link_idx,
                               bool& keep_going,
                               bool& found_cycle) {
        
        auto& curr_link = chain_links[curr_link_idx];
        bool go_to_left = (extending_backward == left_to_right);
        if ((!curr_link.has_left_side && go_to_left) || (!curr_link.has_right_side && !go_to_left)) {
            if (debug) {
                cerr << "reached end of chain" << endl;
            }
            keep_going = false;
            return;
        }
        
        handle_t adjacent_side = go_to_left ? curr_link.left_side : curr_link.right_side;
        if (debug) {
            cerr << "moving to bridge following link " << curr_link_idx << ", leaving out of " << (go_to_left ? "left" : "right") << " side of link at " << graph.get_id(adjacent_side) << (graph.get_is_reverse(adjacent_side) ? "-" : "+") << endl;
            
        }
        
        // add the nodes of the bridge
        size_t bridge_idx = boundary_to_bridge.at(adjacent_side);
        auto& unipath_bridge = unipath_bridges[bridge_idx];
        if (debug) {
            cerr << "adding brige with index " << bridge_idx << " consisting of node(s)" << endl;
            for (auto handle : unipath_bridge) {
                cerr << " " << graph.get_id(handle) << (graph.get_is_reverse(handle) ? "-" : "+") << " (" << id_map.get_name(graph.get_id(handle)) << ")";
            }
            cerr << endl;
        }
        handle_t final_outward; // the last node of the bridge, facing into the next component
        if (unipath_bridge.back() == adjacent_side) {
            // we're traversing the bridge backwards
            if (debug) {
                cerr << "iterating backwards" << endl;
            }
            for (int64_t j = unipath_bridge.size() - 2; j >= 0; --j) {
                if (extending_backward) {
                    graph.prepend_step(curr_paths[0], unipath_bridge[j]);
                    graph.prepend_step(curr_paths[1], unipath_bridge[j]);
                }
                else {
                    graph.append_step(curr_paths[0], graph.flip(unipath_bridge[j]));
                    graph.append_step(curr_paths[1], graph.flip(unipath_bridge[j]));
                }
            }
            final_outward = graph.flip(unipath_bridge.front());
        }
        else {
            // we're traversing the bridge forwards
            assert(adjacent_side == graph.flip(unipath_bridge.front()));
            if (debug) {
                cerr << "iterating forwards" << endl;
            }
            for (size_t j = 1; j < unipath_bridge.size(); ++j) {
                if (extending_backward) {
                    graph.prepend_step(curr_paths[0], graph.flip(unipath_bridge[j]));
                    graph.prepend_step(curr_paths[1], graph.flip(unipath_bridge[j]));
                }
                else {
                    graph.append_step(curr_paths[0], unipath_bridge[j]);
                    graph.append_step(curr_paths[1], unipath_bridge[j]);
                }
            }
            final_outward = unipath_bridge.back();
        }
        if (debug) {
            cerr << "finished adding bridge, current allele paths:" << endl;
            for (auto curr_path : curr_paths) {
                for (auto handle : graph.scan_path(curr_path)) {
                    cerr << " " << graph.get_id(handle) << (graph.get_is_reverse(handle) ? "-" : "+") << " (" << id_map.get_name(graph.get_id(handle)) << ")";
                }
                cerr << endl;
            }
        }
        
        // gather information about the next chain link and its orientation
        auto it = boundary_to_chain_link.find(final_outward);
        if (it == boundary_to_chain_link.end()) {
            // there's no neighboring chainable component
            keep_going = false;
            return;
        }
        curr_link_idx = it->second;
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
        left_to_right = (final_outward == next_link.left_side) != extending_backward;
    };
    
    // function to add the (partial) allele(s) for one haplotype to its chain path
    auto add_next_component = [&](int haplotype,
                                  bool extending_backward,
                                  bool left_to_right, // which direction is "forward" relative to the path
                                  path_handle_t& curr_path,
                                  int64_t curr_link_idx) {
        
        
        auto& next_link = chain_links[curr_link_idx];
        bool enter_left = (extending_backward != left_to_right);
        
        if (debug) {
            cerr << "adding haplotype " << haplotype << " for next link at index " << curr_link_idx << " with alleles" << endl;
            for (bool left : {true, false}) {
                cerr << "from " << (left ? "left" : "right") << endl;
                for (auto allele : (left ? next_link.allele_from_left : next_link.allele_from_right)) {
                    for (auto handle : allele) {
                        cerr << " " << graph.get_id(handle) << (graph.get_is_reverse(handle) ? "-" : "+") << " (" << id_map.get_name(graph.get_id(handle)) << ")";
                    }
                    cerr << endl;
                }
            }
            cerr << "extending backward? " << extending_backward << ", entering from left? " << enter_left << ", link is left-to-right? " << left_to_right << endl;
        }
        
        // handles the logic of adding either the left or right allele onto the path
        auto add_allele = [&](const vector<handle_t>& allele, bool extending) {
            if (enter_left) {
                // we traverse this component left to right
                for (size_t i = extending; i < allele.size(); ++i) {
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
                for (int64_t i = allele.size() - 1 - extending; i >= 0; --i) {
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
            if (debug) {
                cerr << "there is only a one-side allele path for haplotype " << haplotype << endl;
            }
            add_allele(next_link.allele_from_left[haplotype], true);
        }
        else {
            // the walk is broken in the middle
            if (debug) {
                cerr << "there is a broken allele path in this component for haplotype " << haplotype << endl;
            }
            if (enter_left) {
                add_allele(next_link.allele_from_left[haplotype], true);
                // end the haplotype path and start a new one
                curr_path = get_next_path(haplotype);
                add_allele(next_link.allele_from_right[haplotype], false);
            }
            else {
                add_allele(next_link.allele_from_right[haplotype], true);



                // end the haplotype path and start a new one
                curr_path = get_next_path(haplotype);
                add_allele(next_link.allele_from_left[haplotype], false);
            }
        }
        
        if (debug) {
            cerr << "after completing component, current allele path:" << endl;
            for (auto handle : graph.scan_path(curr_path)) {
                cerr << " " << graph.get_id(handle) << (graph.get_is_reverse(handle) ? "-" : "+");
            }
            cerr << endl;
        }
        
        is_chained[curr_link_idx] = true;
    };
    
    for (size_t i = 0; i < chain_links.size(); ++i) {
        if (is_chained[i]) {
            continue;
        }
        
        auto& init_link = chain_links[i];
        if (debug) {
            cerr << "starting a chain at link " << i << " with ";
            if (!init_link.has_left_side) {
                cerr << " no boundary nodes" << endl;
            }
            else if (init_link.has_right_side) {
                cerr << "boundaries " << graph.get_id(init_link.left_side) << (graph.get_is_reverse(init_link.left_side) ? "-" : "+") << " and " << graph.get_id(init_link.right_side) << (graph.get_is_reverse(init_link.right_side) ? "-" : "+") << " (" << id_map.get_name(graph.get_id(init_link.left_side)) << " and " << id_map.get_name(graph.get_id(init_link.right_side)) << ")" << endl;
            }
            else {
                cerr << "single boundary " << graph.get_id(init_link.left_side) << (graph.get_is_reverse(init_link.left_side) ? "-" : "+") << " (" << id_map.get_name(graph.get_id(init_link.left_side)) << ")" << endl;
            }
        }
        
        // make paths to extend into haplotypes
        array<path_handle_t, 2> init_paths;
        for (int hap : {0, 1}) {
            init_paths[hap] = get_next_path(hap);
        }
        array<path_handle_t, 2> curr_paths = init_paths;
        for (int hap : {0, 1}) {
            for (handle_t step : init_link.allele_from_left[hap]) {
                graph.append_step(curr_paths[hap], step);
            }
        }
        is_chained[i] = true;
        
        if (debug) {
            cerr << "initial haplotype paths:" << endl;
            for (auto path : curr_paths) {
                for (auto handle : graph.scan_path(path)) {
                    cerr << " " << graph.get_id(handle) << (graph.get_is_reverse(handle) ? "-" : "+") << " (" << id_map.get_name(graph.get_id(handle)) << ")";
                }
                cerr << endl;
            }
            cerr << "extending to the left" << endl;
        }
        
        bool left_to_right = true;
        int64_t link_index = i;
        bool keep_going = init_link.has_left_side;
        bool found_cycle = false;
        
        // start by extending to the left
        while (true) {
            add_next_bridge(true, left_to_right, curr_paths, link_index,
                            keep_going, found_cycle);
            if (!keep_going) {
                break;
            }
            for (int hap : {0, 1}) {
                add_next_component(hap, true, left_to_right, curr_paths[hap], link_index);
            }
            is_chained[link_index] = true;
        }
        
        if (found_cycle) {
            // we looped back around to the same link, but we might not have included its
            // right side
            
            if (debug) {
                cerr << "finishing a chain that consists of a cycle" << endl;
            }
            
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
            
            if (debug) {
                cerr << "extending to the right" << endl;
            }
            
            // re-initialize the rightward traversal
            keep_going = true;
            left_to_right = true;
            link_index = i;
            for (int hap : {0, 1}) {
                if (init_link.allele_from_right[hap].empty()) {
                    // we can continue the initial allele
                    curr_paths[hap] = init_paths[hap];
                }
                else {
                    // we have to start a new path
                    curr_paths[hap] = get_next_path(hap);
                    for (handle_t step : init_link.allele_from_right[hap]) {
                        graph.append_step(curr_paths[hap], step);
                    }
                }
            }
            
            // add links to the right
            while (true) {
                add_next_bridge(false, left_to_right, curr_paths, link_index,
                                keep_going, found_cycle);
                if (!keep_going) {
                    break;
                }
                for (int hap : {0, 1}) {
                    add_next_component(hap, false, left_to_right, curr_paths[hap], link_index);
                }
                is_chained[link_index] = true;
            }
        }
    }
    
    /*
     * Part 4: finish off with some heuristic post-processing to improve contiguity and correctness
     */
    
    // extend into unphaseable components if we can
    extend_unambiguous_phase_paths(graph, contact_graph);
    
    // self loops are ambiguous, so we don't allow them in phased paths
    break_self_looping_phase_paths(graph);
    
    // we might have broken paths enough that some don't actually have any haploid, phased
    // sequence in them, in which case we remove them
    purge_null_phase_paths(graph, contact_graph);
}

string HamiltonianChainer::phase_path_name(int haplotype, int path_id) {
    return "gfase_hap_" +  to_string(haplotype) + "_" + to_string(path_id);
}

pair<vector<handle_t>, vector<handle_t>>
HamiltonianChainer::generate_allelic_semiwalks(const HandleGraph& graph,
                                               const IncrementalIdMap<string>& id_map,
                                               const unordered_set<nid_t>& in_phase_nodes,
                                               const unordered_set<nid_t>& out_phase_nodes,
                                               const unordered_set<handle_t>& starts,
                                               const unordered_set<handle_t>& ends,
                                               bool& resolved_hamiltonian) const {
    
    if (debug) {
        cerr << "finding alleles for in-phase nodes:" << endl;
        for (auto nid : in_phase_nodes) {
            cerr << "\t" << nid << " (" << id_map.get_name(nid) << ")" << endl;
        }
    }
    
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
        
        if (debug) {
            cerr << "full hamiltonian:" << endl;
            for (auto handle : hamiltonian.hamiltonian_path) {
                cerr << " " << graph.get_id(handle) << (graph.get_is_reverse(handle) ? "-" : "+");
            }
            cerr << endl;
            cerr << "unique prefix:" << endl;
            for (auto handle : hamiltonian.unique_prefix) {
                cerr << " " << graph.get_id(handle) << (graph.get_is_reverse(handle) ? "-" : "+");
            }
            cerr << endl;
        }
        
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
                if (debug) {
                    cerr << "reverse unique prefix:" << endl;
                    for (auto handle : rev_hamiltonian.unique_prefix) {
                        cerr << " " << graph.get_id(handle) << (graph.get_is_reverse(handle) ? "-" : "+");
                    }
                    cerr << endl;
                }
                
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
        if (debug) {
            cerr << "could not resolve hamiltonian allele" << endl;
        }
        if (starts.size() == 1) {
            if (debug) {
                cerr << "attempting to find forward unambiguous path" << endl;
            }
            return_val.first = max_unambiguous_path(graph, *starts.begin(), in_phase_nodes);
            if (debug) {
                cerr << "forward max unambiguous path:" << endl;
                for (auto handle : return_val.first) {
                    cerr << " " << graph.get_id(handle) << (graph.get_is_reverse(handle) ? "-" : "+");
                }
                cerr << endl;
            }
        }
        if (ends.size() == 1) {
            if (debug) {
                cerr << "attempting to find reverse unambiguous path" << endl;
            }
            auto rev_suffix = max_unambiguous_path(graph, *rev_starts.begin(), in_phase_nodes);
            for (auto it = rev_suffix.rbegin(); it != rev_suffix.rend(); ++it) {
                return_val.second.push_back(graph.flip(*it));
            }
            if (debug) {
                cerr << "reverse max unambiguous path:" << endl;
                for (auto handle : return_val.second) {
                    cerr << " " << graph.get_id(handle) << (graph.get_is_reverse(handle) ? "-" : "+");
                }
                cerr << endl;
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
        unambiguous = (unambiguous && num_neighbors == 1);
        if (unambiguous) {
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

void HamiltonianChainer::break_self_looping_phase_paths(MutablePathDeletableHandleGraph& graph) {
    
    for (int hap : {0, 1}) {
        vector<path_handle_t> to_remove, to_add;
        for (auto path_handle : phase_paths[hap]) {
            
            // the steps before where we need to break the path
            unordered_set<step_handle_t> breaks;
            
            unordered_map<handle_t, step_handle_t> previous_steps;
            // scan the path
            for (auto step = graph.path_begin(path_handle), end = graph.path_end(path_handle); step != end; step = graph.get_next_step(step)) {
                
                previous_steps[graph.get_handle_of_step(step)] = step;
                // check if any edges backtrack to an earlier step of the path
                graph.follow_edges(graph.get_handle_of_step(step), false, [&](const handle_t& next) {
                    if (previous_steps.count(next)) {
                        if (step != graph.path_back(path_handle)) {
                            breaks.insert(step);
                        }
                        auto prev_step = previous_steps[next];
                        if (prev_step != graph.path_begin(path_handle)) {
                            breaks.insert(graph.get_previous_step(prev_step));
                        }
                        return false;
                    }
                    return true;
                });
            }
            if (!breaks.empty()) {
                // we will delete the old path
                to_remove.push_back(path_handle);
                
                // scan the path again, breaking it up as we go
                path_handle_t curr_path = graph.create_path_handle(phase_path_name(hap, next_path_ids[hap]++));
                to_add.push_back(curr_path);
                for (auto step = graph.path_begin(path_handle), end = graph.path_end(path_handle); step != end; step = graph.get_next_step(step)) {
                    graph.append_step(curr_path, graph.get_handle_of_step(step));
                    if (breaks.count(step)) {
                        // break the path after this step
                        curr_path = graph.create_path_handle(phase_path_name(hap, next_path_ids[hap]++));
                        to_add.push_back(curr_path);
                    }
                }
            }
        }
        for (auto path_handle : to_remove) {
            graph.destroy_path(path_handle);
            phase_paths[hap].erase(path_handle);
        }
        for (auto path_handle : to_add) {
            phase_paths[hap].insert(path_handle);
        }
    }
}

void HamiltonianChainer::extend_unambiguous_phase_paths(MutablePathDeletableHandleGraph& graph,
                                                        const MultiContactGraph& contact_graph) {
    
    for (int hap : {0, 1}) {
        const auto& hap_phase_paths = phase_paths[hap];
        for (auto path_handle : phase_paths[hap]) {
            for (bool beginning : {true, false}) {
                if (debug) {
                    cerr << "attempting to unambiguously extend phase path " << graph.get_path_name(path_handle) << " at the " << (beginning ? "beginning" : "end") << endl;
                }
                
                auto step = beginning ? graph.path_begin(path_handle) : graph.path_back(path_handle);
                if (step == graph.path_end(path_handle) || step == graph.path_front_end(path_handle)) {
                    if (debug) {
                        cerr << "path is empty, skipping" << endl;
                    }
                    continue;
                }
                
                // are there any other phase paths that end here?
                int num_overlapping = 0;
                graph.for_each_step_on_handle(graph.get_handle_of_step(step), [&](const step_handle_t& overlapping) {
                    num_overlapping += hap_phase_paths.count(graph.get_path_handle_of_step(overlapping));
                });
                if (num_overlapping != 1) {
                    // it can't be unambiguous if we don't know which one to extend
                    continue;
                }
                
                while (true) {
                    handle_t to_add = as_handle(-1);
                    bool unambiguous = graph.follow_edges(graph.get_handle_of_step(step), beginning,
                                                          [&](const handle_t& neighbor) {
                        auto node_id = graph.get_id(neighbor);
                        if (!contact_graph.has_node(node_id) || !contact_graph.has_alt(node_id)) {
                            // it can't be an unambiguous walk because there are allowed, non-phased nodes adjacent
                            return false;
                        }
                        int neighbor_hap = contact_graph.get_partition(node_id) == 1 ? 1 : 0;
                        if (neighbor_hap != hap) {
                            // we can ignore nodes that are phased to the other haplotype
                            return true;
                        }
                        bool off_phase_path = graph.for_each_step_on_handle(neighbor,
                                                                            [&](const step_handle_t& neighbor_step) {
                            return !hap_phase_paths.count(graph.get_path_handle_of_step(neighbor_step));
                        });
                        if (!off_phase_path) {
                            // we can't distinguish between extending this phase path and joining it with the adjacent one
                            return false;
                        }
                        if (to_add == as_handle(-1)) {
                            // this is the first in-phase node we've seen
                            to_add = neighbor;
                            return true;
                        }
                        // we've already seen a different in-phase node
                        return false;
                    });
                    if (unambiguous && to_add != as_handle(-1)) {
                        // we found a single, unambiguous extension into this component
                        if (debug) {
                            cerr << "adding unambiguous extension to " << graph.get_id(to_add) << (graph.get_is_reverse(to_add) ? "-" : "+") << endl;
                        }
                        if (beginning) {
                            graph.prepend_step(path_handle, to_add);
                            step = graph.path_begin(path_handle);
                        }
                        else {
                            graph.append_step(path_handle, to_add);
                            step = graph.path_back(path_handle);
                        }
                    }
                    else {
                        // we can't add any more
                        break;
                    }
                }
            }
        }
    }
}

void HamiltonianChainer::purge_null_phase_paths(MutablePathDeletableHandleGraph& graph,
                                                const MultiContactGraph& contact_graph) {
    
    for (int hap : {0, 1}) {
        vector<path_handle_t> to_remove;
        for (auto path_handle : phase_paths[hap]) {
            // look for any step that is a haploid allele
            bool found_phased_step = false;
            for (auto step = graph.path_begin(path_handle), end = graph.path_end(path_handle); step != end && !found_phased_step; step = graph.get_next_step(step)) {
                auto node_id = graph.get_id(graph.get_handle_of_step(step));
                if (contact_graph.has_node(node_id)) {
                    found_phased_step = contact_graph.has_alt(node_id);
                }
            }
            if (!found_phased_step) {
                if (debug) {
                    cerr << "path " << graph.get_path_name(path_handle) << " did not include any phased sequence, erasing" << endl;
                }
                graph.destroy_path(path_handle);
                to_remove.push_back(path_handle);
            }
        }
        for (auto path_handle : to_remove) {
            phase_paths[hap].erase(path_handle);
        }
    }
}

void HamiltonianChainer::write_chaining_results_to_bandage_csv(path output_dir, const IncrementalIdMap<string>& id_map,
                                                               const MultiContactGraph& contact_graph) const {
    
    path output_path = output_dir / "chains.csv";
    ofstream file(output_path);
    
    if (!file.is_open() || !file.good()) {
        throw std::runtime_error("ERROR: could not write to file: " + output_path.string());
    }
    
    const string cap_0 = "Dark Blue";
    const string cap_1 = "Dark Red";
    const string cap_both = "Indigo";
    const string middle_0 = "Cornflower Blue";
    const string middle_1 = "Crimson";
    const string middle_both = "Dark Orchid";
    const string unchained_0 = "Powder Blue";
    const string unchained_1 = "Light Coral";
    
    file << "Name,Color\n";
    path_graph->for_each_handle([&](const handle_t& handle) {
        bool on_hap_0 = false, on_hap_1 = false, on_hap_0_path = false, on_hap_1_path = false, is_cap = false;
        path_graph->for_each_step_on_handle(handle, [&](const step_handle_t& step) {
            auto path = path_graph->get_path_handle_of_step(step);
            if (phase_paths[0].count(path)) {
                on_hap_0_path = true;
                on_hap_0 = true;
                if (step == path_graph->path_begin(path) || step == path_graph->path_back(path)) {
                    is_cap = true;
                }
            }
            else if (phase_paths[1].count(path)) {
                on_hap_1_path = true;
                on_hap_1 = true;
                if (step == path_graph->path_begin(path) || step == path_graph->path_back(path)) {
                    is_cap = true;
                }
            }
        });
        if (contact_graph.has_node(path_graph->get_id(handle))) {
            if (contact_graph.get_partition(path_graph->get_id(handle)) == -1) {
                on_hap_0 = true;
            }
            else if (contact_graph.get_partition(path_graph->get_id(handle)) == 1) {
                on_hap_1 = true;
            }
        }
        if (debug) {
            cerr << path_graph->get_id(handle) << " / " << id_map.get_name(path_graph->get_id(handle)) << " is on hap0? " << on_hap_0 << ", is on hap1? " << on_hap_1 << " is on hap0 path? " << on_hap_0_path << ", is on hap1 path? " << on_hap_1_path << ", is a cap? " << is_cap << endl;;
        }
        if (on_hap_0 || on_hap_1) {
            string color;
            if (on_hap_0_path || on_hap_1_path) {
                if (is_cap) {
                    if (on_hap_0_path && on_hap_1_path) {
                        color = cap_both;
                    }
                    else if (on_hap_0_path) {
                        color = cap_0;
                    }
                    else {
                        color = cap_1;
                    }
                }
                else {
                    if (on_hap_0_path && on_hap_1_path) {
                        color = middle_both;
                    }
                    else if (on_hap_0_path) {
                        color = middle_0;
                    }
                    else {
                        color = middle_1;
                    }
                }
            }
            else if (on_hap_0) {
                color = unchained_0;
            }
            else {
                color = unchained_1;
            }
            file << id_map.get_name(path_graph->get_id(handle)) << ',' << color << '\n';
        }
    });
}

}
