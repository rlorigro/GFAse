#include "HamiltonianPath.hpp"

#include <unordered_map>
#include <limits>
#include <bitset>
#include <utility>

using std::unordered_map;
using std::numeric_limits;
using std::function;
using std::pair;
using std::bitset;
using std::make_pair;

namespace gfase {

HamilonianProblemResult find_hamiltonian_path(const HandleGraph& graph,
                                              const unordered_set<nid_t>& target_nodes,
                                              const unordered_set<nid_t>& prohibited_nodes,
                                              const unordered_set<handle_t>& allowed_starts,
                                              const unordered_set<handle_t>& allowed_ends,
                                              size_t max_iters) {
    
    HamilonianProblemResult result;
    
    unordered_map<handle_t, size_t> handle_number;
    graph.for_each_handle([&](const handle_t& h) {
        if (!prohibited_nodes.count(graph.get_id(h))) {
            return;
        }
        handle_number[h] = handle_number.size();
        handle_number[graph.flip(h)] = handle_number.size();
    });
    
    if (handle_number.size() > 64) {
        // we can't fit the allowed handles into our bitset
        return result;
    }
    
    function<uint64_t(const handle_t&)> bitcode = [&handle_number](const handle_t& h) {
        return uint64_t(1) << handle_number[h];
    };
    function<bool(uint64_t)> contains_reversal = [](uint64_t set) {
        // int with alternating bits with 0 in 1s place
        static const uint64_t altern_0 = 0xaaaaaaaaaaaaaaaa;
        // int with alternating bits with 1 in 1s place
        static const uint64_t altern_1 = 0x5555555555555555;
        // the two strands are always in alternating positions in the sets
        // so this is only non-zero if they have both 1s set
        return (set & altern_1) & ((set & altern_0) >> 1);
    };
    
    vector<unordered_set<pair<uint64_t, handle_t>>> dp_table(1);
    if (allowed_starts.empty()) {
        // we don't allow unnecesarily long paths, so we still prohibit starting
        // at an unrequired node
        for (nid_t node_id : target_nodes) {
            for (bool reverse : {false, true}) {
                handle_t handle = graph.get_handle(node_id, reverse);
                dp_table[0].emplace(bitcode(handle), handle);
            }
        }
    }
    else {
        for (handle_t handle : allowed_starts) {
            dp_table[0].emplace(bitcode(handle), handle);
        }
    }
    
    // stop if we've hit the longest possible hamiltonian (because we don't allow
    // traversing both strands of a node) or when there are no hamiltonian paths of
    // lenght n - 1
    size_t iter_num = 0;
    while (dp_table.size() < handle_number.size() / 2 && !dp_table.back().empty()
           && iter_num < max_iters) {
        
        dp_table.emplace_back();
        
        auto& prev_table = dp_table[dp_table.size() - 2];
        auto& new_table = dp_table[dp_table.size() - 1];
        
        for (const auto& entry : prev_table) {
            bool keep_going = graph.follow_edges(entry.second, false, [&](const handle_t& next) {
                
                // TODO: i actually only really need to look for the opposite bitcode,
                // not any arbitrary reversals...
                uint64_t new_set = entry.first | bitcode(next);
                if (new_set != entry.first && // was already in set
                    !contains_reversal(new_set) && // reverse is in set
                    !prohibited_nodes.count(graph.get_id(next))) {
                    // extension is valid
                    new_table.emplace(new_set, next);
                }
                
                // give up if this goes on too long
                ++iter_num;
                return iter_num < max_iters;
            });
            if (!keep_going) {
                break;
            }
        }
    }
    
    if (iter_num < max_iters) {
        
        // TODO: the popcount instruction might be faster if i switch to 32-bit integers, but
        // that would require me to template all of this out...
        
        uint64_t completion_code = 0;
        for (nid_t node_id : target_nodes) {
            for (bool reverse : {true, false}) {
                completion_code |= bitcode(graph.get_handle(node_id, reverse));
            }
        }
        auto is_complete = [&](uint64_t set) {
            // c++20 has a template popcount function that could be faster than this, and
            // there are also non-compiler-dependent algorithms, e.g.
            // https://stackoverflow.com/questions/109023/count-the-number-of-set-bits-in-a-32-bit-integer
            return bitset<64>(set & completion_code).count() == target_nodes.size();
        };
        
        // traceback routine
        
        // all of the DP entrys that are part of a traceback in each step
        vector<unordered_set<pair<uint64_t, handle_t>>> traceback_clouds(dp_table.size());
        // a single traceback
        vector<pair<uint64_t, handle_t>> traceback;
        for (int64_t i = dp_table.size() - 1; i >= 0; --i) {
            
            auto& table = dp_table[i];
            auto& cloud = traceback_clouds[i];
            
            // find complete, valid hamiltonian paths
            // we only allow the path to end in a non-target node if that node was provided as a
            // an allowed end (this prevents unnecessary elongation into non-target nodes)
            // it must also not be a part of a traceback that we're already extending
            for (const auto& entry : table) {
                if ((allowed_ends.empty() || allowed_ends.count(entry.second)) &&
                    (allowed_ends.count(entry.second) || target_nodes.count(graph.get_id(entry.second))) &&
                    !cloud.count(entry) &&
                    is_complete(entry.first)) {
                    
                    traceback_clouds[i].emplace(entry);
                    
                    // we always start the single traceback over if we find a new ending, this ensures
                    // that we get the shortest possible path
                    traceback.clear();
                    traceback.push_back(entry);
                }
            }
            
            if (i > 0) {
                // try to find a predecessor for the ongoing tracebacks
                auto& prev_table = dp_table[i - 1];
                auto& prev_cloud = traceback_clouds[i - 1];
                
                for (const auto& dp_entry : cloud) {
                    // remove the final node's bit
                    uint64_t prev_set = dp_entry.first ^ bitcode(dp_entry.second);
                    graph.follow_edges(dp_entry.second, true, [&](const handle_t& prev) {
                        auto prev_entry = make_pair(prev_set, prev);
                        if (prev_table.count(prev_entry)) {
                            prev_cloud.emplace(prev_entry);
                        }
                        if (!traceback.empty() && dp_entry == traceback.back()) {
                            traceback.push_back(prev_entry);
                        }
                    });
                }
            }
        }
        
        // FIXME: we will never get a unique hamiltonian if we don't specify any starts or ends
        // because the orientation can be reversed. we need to break symmetry by insisting on
        // a particular strand for one required node
        
        result.is_solved = true;
        if (!traceback.empty()) {
            // populate the path in the result
            bool found_non_unique = false;
            size_t i = traceback.size() - 1;
            for (auto it = traceback.rbegin(); it != traceback.rend(); ++it) {
                
                result.hamiltonian_path.push_back(it->second);
                
                // check for uniqueness
                found_non_unique = found_non_unique || traceback_clouds[i].size() <= 1;
                if (!found_non_unique) {
                    result.unique_prefix.push_back(it->second);
                }
                --i;
            }
        }
        
    }
    
    // we skip to the end if we run into the maximum number of iterations in the inner loop
    return result;
}

}
