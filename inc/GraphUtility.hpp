#ifndef GFASE_GRAPHUTILITY_HPP
#define GFASE_GRAPHUTILITY_HPP

#include "IncrementalIdMap.hpp"
#include "gfa_to_handle.hpp"
#include "handle_to_gfa.hpp"
#include "Filesystem.hpp"
#include "GfaReader.hpp"
#include "CLI11.hpp"

#include "bdsg/hash_graph.hpp"

#include <functional>
#include <string>
#include <queue>

using gfase::IncrementalIdMap;
using gfase::handle_graph_to_gfa;

using ghc::filesystem::path;
using bdsg::MutablePathMutableHandleGraph;
using bdsg::path_handle_t;
using bdsg::step_handle_t;
using bdsg::handle_t;
using bdsg::HashGraph;

using std::string;
using std::queue;
using std::cout;
using std::cerr;

namespace gfase{

void for_node_in_bfs(HandleGraph& graph, nid_t start_node, const function<void(const handle_t&)>& f);

void for_edge_in_bfs(HandleGraph& graph, nid_t start_node, const function<void(const handle_t& handle_a, const handle_t& handle_b)>& f);

void for_each_connected_component(HandleGraph& graph, const function<void(unordered_set<nid_t>& connected_component)>& f);

}

#endif //GFASE_GRAPHUTILITY_HPP
