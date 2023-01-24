#include "HamiltonianPath.hpp"
#include "IncrementalIdMap.hpp"
#include "gfa_to_handle.hpp"
#include "handle_to_gfa.hpp"
#include "graph_utility.hpp"
#include "Filesystem.hpp"

#include "bdsg/hash_graph.hpp"

#include <string>
#include <set>
#include <iostream>
#include <unordered_set>
#include <algorithm>
#include <utility>
#include <limits>
#include <sstream>

using namespace gfase;

using ghc::filesystem::path;
using bdsg::HashGraph;
using bdsg::HandleGraph;
using handlegraph::handle_t;
using handlegraph::nid_t;

using std::pair;
using std::sort;
using std::string;
using std::cerr;
using std::endl;
using std::unordered_set;
using std::numeric_limits;
using std::stringstream;

string path_to_string(const vector<handle_t>& path, const HandleGraph& graph,
                      const IncrementalIdMap<string>& id_map) {
    stringstream strm;
    for (auto h : path) {
        strm << "\t" << id_map.get_name(graph.get_id(h)) << (graph.get_is_reverse(h) ? "-" : "+") << endl;
    }
    return strm.str();
}

void run_test(string& data_file,
              const unordered_set<pair<string, bool>>& starts,
              const unordered_set<pair<string, bool>>& ends,
              const unordered_set<string>& targets,
              const unordered_set<string>& prohibited,
              bool should_complete,
              size_t minimum_hamiltonian_length,
              size_t unique_prefix_len,
              size_t max_iters = numeric_limits<size_t>::max()) {
    
    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();
    path relative_gfa_path = data_file;
    path absolute_gfa_path = project_directory / relative_gfa_path;
    
    GfaReader reader(absolute_gfa_path);
    
    HashGraph graph;
    IncrementalIdMap<string> id_map;
    
    gfa_to_handle_graph(graph, id_map, absolute_gfa_path);
    
    unordered_set<nid_t> target_nodes, prohibited_nodes;
    for (auto target : targets) {
        target_nodes.insert(id_map.get_id(target));
    }
    for (auto antitarget : prohibited) {
        prohibited_nodes.insert(id_map.get_id(antitarget));
    }
    unordered_set<handle_t> allowed_starts, allowed_ends;
    for (auto start : starts) {
        allowed_starts.insert(graph.get_handle(id_map.get_id(start.first), start.second));
    }
    for (auto end : ends) {
        allowed_ends.insert(graph.get_handle(id_map.get_id(end.first), end.second));
    }
    
    auto result = find_hamiltonian_path(graph, target_nodes, prohibited_nodes,
                                        allowed_starts, allowed_ends, max_iters);
    
    if (!should_complete) {
        if (result.is_solved) {
            stringstream strm;
            strm << "ERROR: Hamiltonian path algorithm exited successfully when it should fail on GFA " << data_file << endl;
            strm << "Hamiltonian path identified:" << endl;
            strm << path_to_string(result.hamiltonian_path, graph, id_map);
            throw runtime_error(strm.str());
        }
    }
    else {
        if (!result.is_solved) {
            stringstream strm;
            strm << "ERROR: Hamiltonian path algorithm failed when it should fail on GFA " << data_file << endl;
            throw runtime_error(strm.str());
        }
        else {
            if (!targets.empty()) {
                const auto& path = result.hamiltonian_path;
                const auto& prefix = result.unique_prefix;
                
                if (path.empty()) {
                    stringstream strm;
                    strm << "ERROR: Empty path found despite specifying target nodes on GFA " << data_file << endl;
                    throw runtime_error(strm.str());
                }
                
                if (!starts.empty()) {
                    pair<string, bool> start(id_map.get_name(graph.get_id(path.front())),
                                             graph.get_is_reverse(path.front()));
                    if (!starts.count(start)) {
                        stringstream strm;
                        strm << "ERROR: Hamiltonian path algorithm does not start at an allowed start on GFA " << data_file << endl;
                        strm << "Hamiltonian path identified:" << endl;
                        strm << path_to_string(path, graph, id_map);
                        throw runtime_error(strm.str());
                    }
                }
                if (!ends.empty()) {
                    pair<string, bool> end(id_map.get_name(graph.get_id(path.back())),
                                           graph.get_is_reverse(path.back()));
                    if (!ends.count(end)) {
                        stringstream strm;
                        strm << "ERROR: Hamiltonian path algorithm does not end at an allowed end on GFA " << data_file << endl;
                        strm << "Hamiltonian path identified:" << endl;
                        strm << path_to_string(path, graph, id_map);
                        throw runtime_error(strm.str());
                    }
                }
                
                unordered_set<nid_t> path_ids;
                for (size_t i = 0; i < path.size(); ++i) {
                    auto nid = graph.get_id(path[i]);
                    if (prohibited_nodes.count(nid)) {
                        stringstream strm;
                        strm << "ERROR: Hamiltonian path contains prohibited node " << id_map.get_name(graph.get_id(path[i])) << " on GFA " << data_file << endl;
                        strm << "Hamiltonian path:" << endl;
                        strm << path_to_string(path, graph, id_map);
                        throw runtime_error(strm.str());
                    }
                    if (path_ids.count(nid)) {
                        stringstream strm;
                        strm << "ERROR: Hamiltonian path contains a repeated node " << id_map.get_name(nid) << " on GFA " << data_file << endl;
                        strm << "Hamiltonian path:" << endl;
                        strm << path_to_string(path, graph, id_map);
                        throw runtime_error(strm.str());
                    }
                    path_ids.insert(nid);
                    if (i > 0) {
                        if (!graph.has_edge(path[i - 1], path[i])) {
                            stringstream strm;
                            strm << "ERROR: Hamiltonian path contains an invalid edge between indexes " << (i - 1) << " and " << i << " on GFA " << data_file << endl;
                            strm << "Hamiltonian path:" << endl;
                            strm << path_to_string(path, graph, id_map);
                            throw runtime_error(strm.str());
                        }
                    }
                }
                for (auto target_name : targets) {
                    if (!path_ids.count(id_map.get_id(target_name))) {
                        stringstream strm;
                        strm << "ERROR: Hamiltonian path is missing target node " << target_name << " on GFA " << data_file << endl;
                        strm << "Hamiltonian path:" << endl;
                        strm << path_to_string(path, graph, id_map);
                        throw runtime_error(strm.str());
                    }
                }
                
                if (path.size() != minimum_hamiltonian_length) {
                    stringstream strm;
                    strm << "ERROR: Hamiltonian is length " << path.size() << " when " << minimum_hamiltonian_length << " was expected on GFA " << data_file << endl;
                    strm << "Hamiltonian path:" << endl;
                    strm << path_to_string(path, graph, id_map);
                    throw runtime_error(strm.str());
                }
                
                if (prefix.size() != unique_prefix_len) {
                    stringstream strm;
                    strm << "ERROR: Unique prefix is length " << prefix.size() << " when " << unique_prefix_len << " was expected on GFA " << data_file << endl;
                    strm << "Hamiltonian path:" << endl;
                    strm << path_to_string(path, graph, id_map);
                    strm << "Unique prefix:" << endl;
                    strm << path_to_string(prefix, graph, id_map);
                    throw runtime_error(strm.str());
                }
                for (size_t i = 0; i < prefix.size(); ++i) {
                    if (prefix[i] != path[i]) {
                        stringstream strm;
                        strm << "ERROR: Unique prefix does not match full path on GFA " << data_file << endl;
                        strm << "Hamiltonian path:" << endl;
                        strm << path_to_string(path, graph, id_map);
                        strm << "Unique prefix:" << endl;
                        strm << path_to_string(prefix, graph, id_map);
                        throw runtime_error(strm.str());
                    }
                }
                
                // TODO: checking the uniqueness of the prefix is actually pretty hard...
            }
            else {
                if (!result.hamiltonian_path.empty()) {
                    stringstream strm;
                    strm << "ERROR: Non-empty path found despite specifying no target nodes on GFA " << data_file << endl;
                    throw runtime_error(strm.str());
                }
            }
        }
    }
}

int main(){
    
    {
        string data_file = "data/simple_chain.gfa";
        
        unordered_set<pair<string, bool>> starts{
            pair<string, bool>("i", false)
        };
        unordered_set<pair<string, bool>> ends{
            pair<string, bool>("m", false)
        };
        unordered_set<string> target_nodes{"a", "d", "e", "h"};
        unordered_set<string> prohibited_nodes{};
        bool should_complete = true;
        size_t minimum_hamiltonian_length = 9;
        size_t unique_prefix_len = 9;
        
        run_test(data_file, starts, ends, target_nodes, prohibited_nodes,
                 should_complete, minimum_hamiltonian_length, unique_prefix_len);
    }
    
    return 0;
}
