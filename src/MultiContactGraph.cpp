#include "MultiContactGraph.hpp"
#include "binomial.hpp"

#include <thread>
#include <ostream>
#include <queue>

using std::thread;
using std::priority_queue;
using std::numeric_limits;
using std::exception;
using std::ofstream;
using std::ostream;
using std::set_intersection;
using std::runtime_error;
using std::queue;
using std::cerr;
using std::min;
using std::max;


namespace gfase{


const char* NonBipartiteEdgeException::what() const noexcept{
    return message.c_str();
}


NonBipartiteEdgeException::NonBipartiteEdgeException(const alt_component_t& c_a, const alt_component_t& c_b, int32_t a, int32_t b):
        runtime_error(""),
        a(a),
        b(b)
{
    component_a = c_a;
    component_b = c_b;

    message += "ERROR: adding alt for " + to_string(a) + ',' + to_string(b) + " would result in non-bipartite component\n";
    message += "component_a:\n";
    for (auto& item: component_a.first){
        message += "0 " + to_string(item) + "\n";
    }
    for (auto& item: component_a.second){
        message += "1 " + to_string(item) + "\n";
    }
    message += "component_b:\n";
    for (auto& item: component_b.first){
        message += "0 " + to_string(item) + "\n";
    }
    for (auto& item: component_b.second){
        message += "1 " + to_string(item) + "\n";
    }

    // a0 AND b0 == {}
    set_intersection(
            component_a.first.begin(),
            component_a.first.end(),
            component_b.first.begin(),
            component_b.first.end(),
            std::inserter(conflicts_0, conflicts_0.begin()));

    // a1 AND b1 == {}
    set_intersection(
            component_a.second.begin(),
            component_a.second.end(),
            component_b.second.begin(),
            component_b.second.end(),
            std::inserter(conflicts_1, conflicts_1.begin()));

    message += "Conflicts found in 0\n";
    for (auto& item: conflicts_0){
        message += (to_string(item) + '\n');
    }

    message += "Conflicts found in 1\n";
    for (auto& item: conflicts_0){
        message += (to_string(item) + '\n');
    }

}


MultiNode::MultiNode(int8_t partition):
        neighbors(),
        coverage(0),
        length(0),
        partition(partition)
{}


bool MultiNode::has_alt() const{
    return not alts.empty();
}


void MultiContactGraph::for_each_alt(int32_t id, const function<void(int32_t alt_id, MultiNode& alt_node)>& f){
    auto& node = nodes.at(id);

    for (auto& alt_id: node.alts){
        auto& alt = nodes.at(alt_id);

        f(alt_id, alt);
    }
}


void MultiContactGraph::for_each_double_alt(int32_t id, const function<void(int32_t alt_id, MultiNode& alt_node)>& f){
    for_each_alt(id, [&](int32_t alt_id, MultiNode& alt_node){
        for_each_alt(alt_id, [&](int32_t alt_id2, MultiNode& alt_node2){
            f(alt_id2, alt_node2);
        });
    });
}


void MultiContactGraph::assert_component_is_valid(const alt_component_t& component) const{
    vector<int32_t> result;
    set_intersection(
            component.first.begin(),
            component.first.end(),
            component.second.begin(),
            component.second.end(),
            std::inserter(result, result.begin()));

    if (not result.empty()){
        for (auto& item: component.first){
            cerr << 0 << ' ' << item << '\n';
        }
        for (auto& item: component.second){
            cerr << 1 << ' ' << item << '\n';
        }

        throw runtime_error("ERROR: Alt component is non-bipartite. See output above for details.");
    }
}


/// Use BFS on node alts to get connected component that represents a bubble
/// \param id
/// \param validate
/// \param component
void MultiContactGraph::get_alt_component(int32_t id, bool validate, alt_component_t& component) const{
    component = {};

    queue <pair <int32_t,int32_t> > q;
    q.emplace(id,0);

    while(not q.empty()){
        auto& [current_id,distance] = q.front();
        q.pop();

        auto result = nodes.find(current_id);

        if (result == nodes.end()){
            throw runtime_error("ERROR: MultiContactGraph::get_alt_component: nonexistent id while iterating: " + to_string(id));
        }

        auto& node = result->second;

        if (distance % 2 == 0){
            component.first.emplace(current_id);
        }
        else{
            component.second.emplace(current_id);
        }

        for (auto& alt_id: node.alts) {
            auto a_result = component.first.find(alt_id);
            auto b_result = component.second.find(alt_id);

            bool a_found = (a_result != component.first.end());
            bool b_found = (b_result != component.second.end());

            if (not (a_found or b_found)){
                q.emplace(alt_id, distance+1);
            }
        }
    }

    if (validate){
        assert_component_is_valid(component);
    }
}


/// Use BFS on node alts to get connected component that represents a bubble
/// \param id
/// \param validate
/// \param component
bool MultiContactGraph::of_same_component_side(int32_t id_a, int32_t id_b) const{
    set <int32_t> visited;
    queue <pair <int32_t,int32_t> > q;
    q.emplace(id_a,0);

    bool same_side = false;

    while(not q.empty()){
        auto& [current_id,distance] = q.front();
        q.pop();

        auto result = nodes.find(current_id);

        if (result == nodes.end()){
            throw runtime_error("ERROR: MultiContactGraph::of_same_component_side: nonexistent id while iterating: " + to_string(current_id));
        }

        auto& node = result->second;

        if (distance % 2 == 0 and current_id == id_b){
            same_side = true;
            break;
        }

        visited.emplace(current_id);

        for (auto& alt_id: node.alts) {
            auto r = visited.find(alt_id);
            bool found = (r != visited.end());

            if (not found){
                q.emplace(alt_id, distance+1);
            }
        }
    }

    return same_side;
}


/// Use BFS on node alts to get connected component that represents a bubble
/// \param id
/// \param validate
/// \param component
bool MultiContactGraph::of_same_component(int32_t id_a, int32_t id_b) const{
    set <int32_t> visited;
    queue <pair <int32_t,int32_t> > q;
    q.emplace(id_a,0);

    bool same_side = false;

    while(not q.empty()){
        auto& [current_id,distance] = q.front();
        q.pop();

        auto result = nodes.find(current_id);

        if (result == nodes.end()){
            throw runtime_error("ERROR: MultiContactGraph::of_same_component_side: nonexistent id while iterating: " + to_string(current_id));
        }

        auto& node = result->second;

        if (current_id == id_b){
            same_side = true;
            break;
        }

        visited.emplace(current_id);

        for (auto& alt_id: node.alts) {
            auto r = visited.find(alt_id);
            bool found = (r != visited.end());

            if (not found){
                q.emplace(alt_id, distance+1);
            }
        }
    }

    return same_side;
}


const array<string,3> MultiContactGraph::colors = {"Cornflower Blue", "Plum", "Tomato"};


MultiContactGraph::MultiContactGraph(const contact_map_t& contact_map, const IncrementalIdMap<string>& id_map){
//    nodes.set_deleted_key(-1);

    string s;

    for (const auto& [a,sub_map]: contact_map){
        try_insert_node(a);

        for (const auto& [b,count]: sub_map){
            try_insert_node(b);
            insert_edge(a,b,count);
        }
    }
}


MultiContactGraph::MultiContactGraph(path csv_path, const IncrementalIdMap<string>& id_map){
    ifstream file(csv_path);

    if (not (file.is_open() and file.good())){
        throw runtime_error("ERROR: could not read file: " + csv_path.string());
    }

    char c;

    string a;
    string b;
    string weight_token;
    int32_t weight;

    size_t n_delimiters = 0;
    size_t n_lines = 0;

    while (file.get(c)){
        if (c == '\n'){
            if (n_lines > 0){
                auto id_a = int32_t(id_map.get_id(a));
                auto id_b = int32_t(id_map.get_id(b));

                weight = stoi(weight_token);

                try_insert_node(id_a);
                try_insert_node(id_b);
                try_insert_edge(id_a, id_b, weight);
            }

            a.clear();
            b.clear();
            weight_token.clear();

            n_delimiters = 0;
            n_lines++;
        }
        else if (c == ','){
            n_delimiters++;
        }
        else{
            if (n_delimiters == 0){
                a += c;
            }
            else if (n_delimiters == 1){
                b += c;
            }
            else if (n_delimiters == 2){
                weight_token += c;
            }
            else{
                throw runtime_error("ERROR: too many delimiters for line in file: " + csv_path.string());
            }
        }
    }
}


//MultiContactGraph::MultiContactGraph(){
//    nodes.set_deleted_key(-1);
//}


void MultiContactGraph::insert_edge(int32_t a, int32_t b, int32_t weight){
    edge_weights.emplace(edge(a,b),weight);
    nodes.at(a).neighbors.emplace(b);
    nodes.at(b).neighbors.emplace(a);
}


void MultiContactGraph::try_insert_edge(int32_t a, int32_t b){
    if (nodes.count(a) == 0 or nodes.count(b) == 0){
        throw runtime_error("ERROR: cannot add edge with nonexistent node id: (" + to_string(a) + "," + to_string(b) +")");
    }

    auto e = edge(a,b);
    auto result = edge_weights.find(e);

    if (result == edge_weights.end()) {
        insert_edge(a,b,0);
    }
}


bool MultiContactGraph::has_edge(int32_t a, int32_t b) const{
    if (nodes.count(a) == 0 or nodes.count(b) == 0){
        return false;
    }

    bool has_edge = false;

    auto e = edge(a,b);
    auto result = edge_weights.find(e);

    if (result != edge_weights.end()) {
        has_edge = true;
    }

    return has_edge;
}


bool MultiContactGraph::components_are_compatible(const alt_component_t& component_a, const alt_component_t& component_b) const{
    // a0 AND b0 == {}
    vector<int32_t> a0_u_b0;
    set_intersection(
            component_a.first.begin(),
            component_a.first.end(),
            component_b.first.begin(),
            component_b.first.end(),
            std::inserter(a0_u_b0, a0_u_b0.begin()));

    // a1 AND b1 == {}
    vector<int32_t> a1_u_b1;
    set_intersection(
            component_a.second.begin(),
            component_a.second.end(),
            component_b.second.begin(),
            component_b.second.end(),
            std::inserter(a1_u_b1, a1_u_b1.begin()));

    return a0_u_b0.empty() and a1_u_b1.empty();
}


void MultiContactGraph::merge_components(
        const alt_component_t& component_a,
        const alt_component_t& component_b,
        alt_component_t& merged_component
        ) const{

    merged_component.first.clear();
    merged_component.second.clear();

    set_union(
            component_a.first.begin(),
            component_a.first.end(),
            component_b.second.begin(),
            component_b.second.end(),
            std::inserter(merged_component.first, merged_component.first.begin()));

    set_union(
            component_a.second.begin(),
            component_a.second.end(),
            component_b.first.begin(),
            component_b.first.end(),
            std::inserter(merged_component.second, merged_component.second.begin()));
}


void MultiContactGraph::add_alt(const alt_component_t& a, const alt_component_t& b, bool remove_weights) {
    if (components_are_compatible(a, b)){
        alt_component_t merged_component;
        merge_components(a, b, merged_component);

        vector<int32_t> id_list;

        id_list.reserve(merged_component.first.size() + merged_component.second.size());

        for (auto id_a: merged_component.first) {
            id_list.emplace_back(id_a);
        }
        for (auto id_b: merged_component.second) {
            id_list.emplace_back(id_b);
        }

        // No valid weights can exist between nodes of a component
        if (remove_weights){
            for (size_t i=0; i<id_list.size(); i++){
                for (size_t j=i+1; j<id_list.size(); j++){
                    auto id_a = id_list.at(i);
                    auto id_b = id_list.at(j);
                    remove_edge(id_a, id_b);
                }
            }
        }

        for (auto id_a: merged_component.first) {
            auto& node_a = nodes.at(id_a);

            for (auto id_b: merged_component.second) {
                auto& node_b = nodes.at(id_b);

                // Enforce all-vs-all connectivity in alt components
                node_a.alts.emplace(id_b);
                node_b.alts.emplace(id_a);

                // Assign partition
                node_a.partition = 1;
                node_b.partition = -1;
            }
        }
    }
    else{
        NonBipartiteEdgeException e(a, b, -1, -1);
        throw e;
    }
}


void MultiContactGraph::add_alt(int32_t a, int32_t b){
    if (nodes.count(a) == 0 or nodes.count(b) == 0){
        throw runtime_error("ERROR: cannot add alt with nonexistent node id: (" + to_string(a) + "," + to_string(b) +")");
    }

    if (a == b){
        throw runtime_error("ERROR: cannot add alt to itself: " + to_string(b));
    }

    if (nodes.at(a).alts.count(b) > 0){
        return;
    }

    // Start by doing alt-wise BFS to get bipartite component of nodes
    alt_component_t component_a;
    alt_component_t component_b;
    get_alt_component(a, false, component_a);
    get_alt_component(b, false, component_b);

    add_alt(component_a, component_b);
}


void MultiContactGraph::try_insert_edge(int32_t a, int32_t b, int32_t weight){
    if (nodes.count(a) == 0 or nodes.count(b) == 0){
        throw runtime_error("ERROR: cannot add edge with nonexistent node id: (" + to_string(a) + "," + to_string(b) +")");
    }

    auto e = edge(a,b);
    auto result = edge_weights.find(e);

    if (result == edge_weights.end()) {
        insert_edge(a,b,weight);
    }
}


int32_t MultiContactGraph::get_edge_weight(int32_t a, int32_t b) const{
    int32_t weight = 0;

    auto e = edge(a,b);
    auto result = edge_weights.find(e);

    if (result != edge_weights.end()) {
        weight = result->second;
    }

    return weight;
}


void MultiContactGraph::increment_edge_weight(int32_t a, int32_t b, int32_t value){
    if (nodes.count(a) == 0 or nodes.count(b) == 0){
        throw runtime_error("ERROR: cannot add edge with nonexistent node id: (" + to_string(a) + "," + to_string(b) +")");
    }

    auto e = edge(a,b);
    auto result = edge_weights.find(e);

    if (result != edge_weights.end()) {
        result->second += value;
    }
}


void MultiContactGraph::increment_coverage(int32_t id, int64_t value){
    if (nodes.count(id) == 0){
        throw runtime_error("ERROR: cannot update coverage for nonexistent node id: " + to_string(id));
    }

    nodes.at(id).coverage += value;
}


void MultiContactGraph::set_node_coverage(int32_t id, int64_t value){
    if (nodes.count(id) == 0){
        throw runtime_error("ERROR: cannot update coverage for nonexistent node id: " + to_string(id));
    }

    nodes.at(id).coverage = value;
}


void MultiContactGraph::set_node_length(int32_t id, int32_t length){
    if (nodes.count(id) == 0){
        throw runtime_error("ERROR: cannot update length for nonexistent node id: " + to_string(id));
    }

    nodes.at(id).length = length;
}


void MultiContactGraph::remove_edge(int32_t a, int32_t b){
    auto e = edge(a,b);
    auto result = edge_weights.find(e);

    if (result != edge_weights.end()) {
        edge_weights.erase(result);
        nodes.at(a).neighbors.erase(b);
        nodes.at(b).neighbors.erase(a);
    }
}


void MultiContactGraph::for_each_node_neighbor(int32_t id, const function<void(int32_t id_other, const MultiNode& n)>& f) const{
    auto result = nodes.find(id);

    if (result == nodes.end()){
        throw runtime_error("ERROR: MultiContactGraph::for_each_node_neighbor: cannot iterate neighbors for id not in contact graph: " + to_string(id));
    }

    auto& node = result->second;

    for (auto& id_other: node.neighbors){
        auto result_other = nodes.find(id_other);

        if (result_other == nodes.end()){
            throw runtime_error("ERROR: MultiContactGraph::for_each_node_neighbor: cannot find neighbor node for id not in contact graph: " + to_string(id));
        }

        auto& node_other = result_other->second;

        f(id_other, node_other);
    }
}


void MultiContactGraph::for_each_node(const function<void(int32_t id, const MultiNode& n)>& f) const{
    for (auto& [id,node]: nodes){
        f(id, node);
    }
}


void MultiContactGraph::for_each_node(const function<void(int32_t id)>& f) const{
    for (auto& [id,node]: nodes){
        f(id);
    }
}


void MultiContactGraph::for_each_edge(const function<void(const pair<int32_t,int32_t> edge, int32_t weight)>& f) const{
    for (auto& [e, weight]: edge_weights){
        f(e, weight);
    }
}


void MultiContactGraph::insert_node(int32_t id, int8_t partition){
    if (partition < -1 or partition > 1){
        throw runtime_error("ERROR: can't assign partition index outside of {-1,0,1}");
    }

    nodes.emplace(id, partition);
}


void MultiContactGraph::insert_node(int32_t id){
    nodes.emplace(id, 0);
}


void MultiContactGraph::try_insert_node(int32_t id){
    if (nodes.count(id) == 0) {
        nodes.emplace(id, 0);
    }
}


void MultiContactGraph::try_insert_node(int32_t id, int8_t partition){
    if (nodes.count(id) == 0) {
        nodes.emplace(id, partition);
    }
}


void MultiContactGraph::validate_alts() {
    for (auto& [id,node]: nodes){

        // If this node is linked to an alt, alt must be maintained in an opposite state
        if (node.has_alt()){
            for (auto& alt_id: node.alts){
                auto& alt = nodes.at(alt_id);

                cerr << id << ',' << alt_id << ',' << int(node.partition) << ',' << int(alt.partition) << '\n';

                if (alt.partition == node.partition){
                    throw runtime_error("ERROR: (MultiContactGraph::set_partition) alt nodes in same partition: " + to_string(int(id)) + ',' + to_string(int(alt_id)));
                }
            }
        }
    }
}


int8_t MultiContactGraph::get_partition(int32_t id) const{
    auto result = nodes.find(id);

    if (result == nodes.end()){
        throw runtime_error("ERROR: cannot find partition for nonexistent node ID: " + to_string(id));
    }

    return result->second.partition;
}


void MultiContactGraph::set_partition(int32_t id, int8_t partition) {
    auto result = nodes.find(id);

    if (result == nodes.end()){
        throw runtime_error("ERROR: MultiContactGraph::set_partition: cannot set partition for id not in contact graph: " + to_string(id));
    }

    auto& node = result->second;
    node.partition = partition;

    // If this node is linked to an alt, alt must be maintained in an opposite state,
    // and double-alts must be maintained in identical state
    if (node.has_alt()){
        if (partition == 0) {
            throw runtime_error("ERROR: cannot set 0 partition for bubble: " + to_string(id));
        }

        alt_component_t component;
        get_alt_component(id, false, component);

        for (auto& alt_id: component.second){
            auto& alt = nodes.at(alt_id);
            alt.partition = partition*int8_t(-1);
        }
        for (auto& alt_id: component.first){
            auto& alt = nodes.at(alt_id);
            alt.partition = partition;
        }
    }
}


void MultiContactGraph::set_partition(const alt_component_t& component, int8_t partition) {
    for (auto& id: component.first){
        auto& a = nodes.at(id);
        a.partition = partition;
    }
    for (auto& id: component.second){
        auto& b = nodes.at(id);
        b.partition = partition*int8_t(-1);
    }
}


void MultiContactGraph::remove_node(int32_t id){
    vector <pair <int32_t, int32_t> > to_be_removed;

    for_each_node_neighbor(id, [&](int32_t id_other, const MultiNode& n){
        to_be_removed.emplace_back(edge(id, id_other));
    });

    for (auto& e: to_be_removed){
        remove_edge(e.first, e.second);
    }

    // Make sure there is no dangling reference to this node in its alt
    auto& n = nodes.at(id);
    if (n.has_alt()) {
        for (auto& alt_id: n.alts){
            nodes.at(alt_id).alts.erase(id);
        }
    }

    nodes.erase(id);
}


size_t MultiContactGraph::edge_count(int32_t id){
    return nodes.at(id).neighbors.size();
}


size_t MultiContactGraph::edge_count(){
    return edge_weights.size();
}


size_t MultiContactGraph::size(){
    return nodes.size();
}


//void MultiContactGraph::resize(){
//    nodes.resize(0);
//}


double MultiContactGraph::get_score(const MultiNode& a, const MultiNode& b, int32_t weight) const{
    double score = 0;

    auto p_a = a.partition;
    auto p_b = b.partition;

    if (p_a != 0 and p_b != 0) {
        score = p_a * p_b * weight;
    }
    else{
        score = 0;
    }

    return score;
}


double MultiContactGraph::compute_consistency_score(int32_t id) const{
    double score = 0;

    auto n = nodes.at(id);

    for_each_node_neighbor(id, [&](int32_t id_other, const MultiNode& n_other){
        // Skip self edges if there are any
        if (id == id_other){
            return;
        }

        score += get_score(n, n_other, edge_weights.at(edge(id, id_other)));
//        cerr << '\t' << id << "<->" << id_other << ' ' << int(n.partition) << 'x' << int(n_other.partition) << 'x' << edge_weights.at(edge(id, id_other)) << ' ' << score << '\n';
    });

    for (auto& alt_id: n.alts){
        auto& n_alt = nodes.at(alt_id);
        for_each_node_neighbor(alt_id, [&](int32_t id_other, const MultiNode& n_other){
            score += get_score(n_alt, n_other, edge_weights.at(edge(alt_id, id_other)));
        });
    }

    return score;
}


double MultiContactGraph::compute_consistency_score(alt_component_t& component) const{
    double score = 0;

    for (auto& id: component.first){
        auto n = nodes.at(id);

        for_each_node_neighbor(id, [&](int32_t id_other, const MultiNode& n_other){
            // Skip self edges if there are any
            if (id == id_other){
                return;
            }

            score += get_score(n, n_other, edge_weights.at(edge(id, id_other)));
//            cerr << '\t' << id << "<->" << id_other << ' ' << int(n.partition) << 'x' << int(n_other.partition) << 'x' << edge_weights.at(edge(id, id_other)) << ' ' << score << '\n';
        });
    }

    for (auto& id: component.second){
        auto n = nodes.at(id);

        for_each_node_neighbor(id, [&](int32_t id_other, const MultiNode& n_other){
            // Skip self edges if there are any
            if (id == id_other){
                return;
            }

            score += get_score(n, n_other, edge_weights.at(edge(id, id_other)));
//            cerr << '\t' << id << "<->" << id_other << ' ' << int(n.partition) << 'x' << int(n_other.partition) << 'x' << edge_weights.at(edge(id, id_other)) << ' ' << score << '\n';
        });
    }

    return score;
}


double MultiContactGraph::compute_total_consistency_score() const{
    double score = 0;

    for_each_edge([&](const pair<int32_t,int32_t> edge, int32_t weight){
        // Skip self edges if any exist
        if (edge.first == edge.second){
            return;
        }

        auto& a = nodes.at(edge.first);
        auto& b = nodes.at(edge.second);
        score += get_score(a, b, weight);
    });

    return score;
}


void MultiContactGraph::randomize_partitions(){
    // True random number
    std::random_device rd;

    // Pseudorandom generator with true random seed
    std::mt19937 rng(rd());
    std::uniform_int_distribution<int> uniform_distribution(0,2);

    for (const auto& [id,node]: nodes){
        int8_t p;
        if (node.has_alt()){
            // Only allow {1,-1} for known bubbles
            p = int8_t((uniform_distribution(rng) % 2));

            if (p == 0){
                p = -1;
            }

            set_partition(id,p);
        }
        else{
            // Allow {1,0,-1}
            p = int8_t((uniform_distribution(rng) % 3) - 1);

            set_partition(id,p);
        }
    }
}


void MultiContactGraph::get_alt_components(vector <alt_component_t>& alt_components) const{
    alt_components.clear();

    unordered_set<int32_t> visited;
    visited.reserve(nodes.size());
    alt_component_t component;

    for (const auto& [n,node]: nodes){
        if (visited.count(n) > 0){
            continue;
        }

        get_alt_component(n, false, component);

        alt_components.emplace_back(component);

        // Dump all the ids from both sides of the component into the list of visited nodes
        for (auto id: component.first){
            visited.insert(id);
        }
        for (auto id: component.second){
            visited.insert(id);
        }
    }
}


void MultiContactGraph::get_alt_component_representatives(vector<int32_t>& representative_ids) const{
    unordered_set<int32_t> visited;
    visited.reserve(nodes.size());
    alt_component_t component;

    for (const auto& [n,node]: nodes){
        if (visited.count(n) > 0){
            continue;
        }

        representative_ids.emplace_back(n);

        get_alt_component(n, false, component);

        // Dump all the ids from both sides of the component into the list of visited nodes
        for (auto id: component.first){
            visited.insert(id);
        }
        for (auto id: component.second){
            visited.insert(id);
        }
    }
}


void MultiContactGraph::get_partitions(vector <pair <int32_t,int8_t> >& partitions) const{
    partitions.clear();
    partitions.reserve(nodes.size());

    for (const auto& [n, node]: nodes){
        partitions.emplace_back(n,node.partition);
    }
}


void MultiContactGraph::set_partitions(const vector <pair <int32_t,int8_t> >& partitions){
    for (const auto& [n, p]: partitions){
        set_partition(n, p);
    }
}


void MultiContactGraph::get_node_ids(vector<int32_t>& ids){
    ids.reserve(nodes.size());

    for (const auto& [id,node]: nodes){
        ids.emplace_back(id);
    }
}


bool MultiContactGraph::has_alt(int32_t id) const{
    auto result = nodes.find(id);

    if (result == nodes.end()){
        throw runtime_error("ERROR: cannot find alt for id not in contact graph: " + to_string(id));
    }

    return result->second.has_alt();
}


bool MultiContactGraph::has_node(int32_t id) const{
    return nodes.find(id) != nodes.end();
}


int64_t MultiContactGraph::get_node_coverage(int32_t id) const{
    auto result = nodes.find(id);

    if (result == nodes.end()){
        throw runtime_error("ERROR: cannot get coverage for nonexistent node: " + to_string(id));
    }

    return result->second.coverage;
}


int32_t MultiContactGraph::get_node_length(int32_t id) const{
    auto result = nodes.find(id);

    if (result == nodes.end()){
        throw runtime_error("ERROR: cannot get coverage for nonexistent node: " + to_string(id));
    }

    return result->second.length;
}


ostream& operator<<(ostream& o, const MultiNode& n){
    o << '\t' << "partition: " << int(n.partition) << '\n';
    o << '\t' << "neighbors: ";

    for (auto& id: n.neighbors){
        o << id << ' ';
    }

    return o;
}


void MultiContactGraph::write_bandage_csv(path output_path, const IncrementalIdMap<string>& id_map) const{
    ofstream file(output_path);

    if (not file.is_open() or not file.good()) {
        throw std::runtime_error("ERROR: could not write to file: " + output_path.string());
    }

    file << "Name" << ',' << "Phase" << ',' << "Coverage" << ',' << "Length" << ',' << "Color" << '\n';

    for (const auto& [id,node]: nodes){
        auto name = id_map.get_name(id);
        file << name << ',' << int(node.partition) << ',' << node.coverage << ',' << node.length << ',' << colors[node.partition+1] << '\n';
    }
}


void MultiContactGraph::write_node_data(path output_path, const IncrementalIdMap<string>& id_map) const{
    ofstream file(output_path);

    if (not file.is_open() or not file.good()) {
        throw std::runtime_error("ERROR: could not write to file: " + output_path.string());
    }

    file << "Id" << ',' << "Name" << ',' << "Coverage" << ',' << "Length" << '\n';

    for (auto& [id,node]: nodes){
        auto name = id_map.get_name(id);
        file << id << ',' << name << ',' << node.coverage << ',' << node.length << '\n';
    }
}


void MultiContactGraph::write_contact_map(path output_path, const IncrementalIdMap<string>& id_map) const{
    ofstream output_file(output_path);

    if (not output_file.is_open() or not output_file.good()) {
        throw std::runtime_error("ERROR: could not write to file: " + output_path.string());
    }

    output_file << "name_a" << ',' << "name_b" << ',' << "weight" << '\n';

    for_each_edge([&](const pair<int32_t, int32_t> edge, int32_t weight) {
        output_file << id_map.get_name(edge.first) << ',' << id_map.get_name(edge.second) << ',' << weight << '\n';
    });
}


void random_node_scale_phase_search(
        MultiContactGraph contact_graph,
        const vector<int32_t>& ids,
        vector <pair <int32_t,int8_t> >& best_partitions,
        atomic<double>& best_score,
        atomic<size_t>& job_index,
        mutex& phase_mutex,
        size_t m_iterations
){

    size_t m = job_index.fetch_add(1);

    // True random number
    std::random_device rd;

    // Pseudorandom generator with true random seed
    std::mt19937 rng(rd());
    std::uniform_int_distribution<int> uniform_distribution(0,int(ids.size()-1));

    contact_graph.set_partitions(best_partitions);

    double total_score;

    while (m < m_iterations) {
        // Randomly perturb
        for (size_t i=0; i<((ids.size()/30) + 1); i++) {
            auto r = ids.at(uniform_distribution(rng));

            int8_t p;
            if (contact_graph.has_alt(r)){
                // Only allow {1,-1}
                p = int8_t((uniform_distribution(rng) % 2));

                if (p == 0){
                    p = -1;
                }
            }
            else{
                // Allow {1,0,-1}
                p = int8_t((uniform_distribution(rng) % 3) - 1);
            }

            contact_graph.set_partition(r, p);
        }

        for (size_t i=0; i<ids.size(); i++) {
            auto n = ids.at(uniform_distribution(rng));

            int64_t max_score = std::numeric_limits<int64_t>::min();
            int8_t p_max = 0;

            if (contact_graph.edge_count(n) == 0){
                continue;
            }

            bool has_alt = contact_graph.has_alt(n);

            for (int8_t p=-1; p<=1; p++) {
                // If the node has an "alt" it can't be made neutral
                if (has_alt and p==0){
                    continue;
                }

                contact_graph.set_partition(n, p);
                auto score = contact_graph.compute_consistency_score(n);

                if (score > max_score) {
                    max_score = score;
                    p_max = p;
                }

//                cerr << n << ' ' << int(has_alt) << ' ' << int(p) << ' ' << score << ' ' << max_score << '\n';
            }

            contact_graph.set_partition(n, p_max);
        }

        total_score = contact_graph.compute_total_consistency_score();

        phase_mutex.lock();
        if (total_score > best_score) {
            best_score = total_score;
            contact_graph.get_partitions(best_partitions);
        }
        else {
            contact_graph.set_partitions(best_partitions);
        }

//        cerr << m << ' ' << best_score << ' ' << total_score << ' ' << std::flush;
//        for (auto& [n,p]: best_partitions){
//            cerr << '(' << n << ',' << int(p) << ") ";
//        }
//        cerr << '\n';

        phase_mutex.unlock();

        m = job_index.fetch_add(1);
    }

    contact_graph.set_partitions(best_partitions);
}



void random_component_scale_phase_search(
        MultiContactGraph contact_graph,
        const vector<int32_t>& ids,
        vector <pair <int32_t,int8_t> >& best_partitions,
        atomic<double>& best_score,
        atomic<size_t>& job_index,
        mutex& phase_mutex,
        size_t m_iterations
){

    size_t m = job_index.fetch_add(1);

    // True random number
    std::random_device rd;

    // Pseudorandom generator with true random seed
    std::mt19937 rng(rd());
    std::uniform_int_distribution<int> uniform_distribution(0,int(ids.size()-1));

    contact_graph.set_partitions(best_partitions);

    alt_component_t component;

    double total_score;

    while (m < m_iterations) {
        // Randomly perturb
        for (size_t i=0; i<((ids.size()/30) + 1); i++) {
            auto r = ids.at(uniform_distribution(rng));

            int8_t p;
            if (contact_graph.has_alt(r)){
                // Only allow {1,-1}
                p = int8_t((uniform_distribution(rng) % 2));

                if (p == 0){
                    p = -1;
                }
            }
            else{
                // Allow {1,0,-1}
                p = int8_t((uniform_distribution(rng) % 3) - 1);
            }

            contact_graph.set_partition(r, p);
        }

        for (size_t i=0; i<ids.size(); i++) {
            auto n = ids.at(uniform_distribution(rng));

            int64_t max_score = std::numeric_limits<int64_t>::min();
            int8_t p_max = 0;

            contact_graph.get_alt_component(n, false, component);
            bool has_alt = not component.first.empty() and not component.second.empty();

            for (int8_t p=-1; p<=1; p++) {
                // If the node has an "alt" it can't be made neutral
                if (has_alt and p==0){
                    continue;
                }

                contact_graph.set_partition(component, p);

                auto score = contact_graph.compute_consistency_score(component);

                if (score > max_score) {
                    max_score = score;
                    p_max = p;
                }

//                cerr << n << ' ' << int(has_alt) << ' ' << int(p) << ' ' << score << ' ' << max_score << '\n';
            }

            contact_graph.set_partition(n, p_max);
        }

        total_score = contact_graph.compute_total_consistency_score();

        phase_mutex.lock();
        if (total_score > best_score) {
            best_score = total_score;
            contact_graph.get_partitions(best_partitions);
        }
        else {
            contact_graph.set_partitions(best_partitions);
        }

//        cerr << m << ' ' << best_score << ' ' << total_score << ' ' << std::flush;
//        for (auto& [n,p]: best_partitions){
//            cerr << '(' << n << ',' << int(p) << ") ";
//        }
//        cerr << '\n';

        phase_mutex.unlock();

        m = job_index.fetch_add(1);
    }

    contact_graph.set_partitions(best_partitions);
}

using orientation_edge_t = pair <int32_t,int32_t>;
using orientation_weight_t = array<int32_t, 2>;

class OrientationDistribution{
public:
    unordered_map <orientation_edge_t, orientation_weight_t> edge_weights;
    vector <alt_component_t> alt_components;

    OrientationDistribution(const MultiContactGraph& contact_graph);
    void write_contact_map(path output_path, const IncrementalIdMap<string>& id_map) const;
    void update(const MultiContactGraph& contact_graph);
};


class OrientationEdgeComparator{
public:
    bool operator()(
            const pair<orientation_edge_t,orientation_weight_t>& a,
            const pair<orientation_edge_t,orientation_weight_t>& b
            ){

        auto a_ordinal = max(a.second[0], a.second[1]);
        auto b_ordinal = max(b.second[0], b.second[1]);

        bool result;

        if (a_ordinal == b_ordinal){
            result = a.second[0] < b.second[0];
        }
        else{
            result = a_ordinal < b_ordinal;
        }

        return result;
    }
};


void OrientationDistribution::write_contact_map(path output_path, const IncrementalIdMap<string>& id_map) const{
    ofstream output_file(output_path);

    if (not output_file.is_open() or not output_file.good()) {
        throw std::runtime_error("ERROR: could not write to file: " + output_path.string());
    }

    output_file << "name_a" << ',' << "name_b" << ',' << "weight_0" << ',' << "weight_1" << ',' << "p_null" << '\n';

    for(auto& [edge, weights]: edge_weights) {
        double n = weights[0] + weights[1];
        double k = min(weights[0], weights[1]);
        output_file << id_map.get_name(edge.first) << ',' << id_map.get_name(edge.second) << ',' << weights[0] << ',' << weights[1] << ',' << binomial(0.5,n,k) << '\n';
    }
}


OrientationDistribution::OrientationDistribution(const MultiContactGraph& contact_graph){
    contact_graph.get_alt_components(alt_components);

    contact_graph.for_each_edge([&](const pair<int32_t,int32_t> edge, int32_t weight){
        edge_weights.insert({edge, {0,0}});
    });
}


void OrientationDistribution::update(const MultiContactGraph& contact_graph){
    contact_graph.for_each_edge([&](const pair<int32_t,int32_t> edge, int32_t weight){
        auto [a,b] = edge;

        bool orientation = contact_graph.get_partition(a) == contact_graph.get_partition(b);

        edge_weights[edge][orientation]++;
    });
}


void flip_component(alt_component_t& c){
    auto temp = c.second;
    c.second = c.first;
    c.first = temp;
}


void monte_carlo_phase_contacts(
        MultiContactGraph& contact_graph,
        const IncrementalIdMap<string>& id_map,
        path output_dir,
        size_t n_threads
){
    auto unmerged_contact_graph = contact_graph;

    size_t m_iterations = 200;
    size_t sample_size = 30;
    size_t n_rounds = 2;

    unordered_set<orientation_edge_t> visited_edges;

    vector <pair <int32_t,int8_t> > phase_state;

    for (size_t i=0; i<n_rounds; i++){
        OrientationDistribution orientation_distribution(contact_graph);

        cerr << "---- " << i << " ----" << '\n';
        for (size_t j=0; j<sample_size; j++){
            phase_contacts(contact_graph, n_threads, m_iterations);
            auto j_score = contact_graph.compute_total_consistency_score();

            contact_graph.get_partitions(phase_state);
            unmerged_contact_graph.set_partitions(phase_state);
            auto j_score_unmerged = unmerged_contact_graph.compute_total_consistency_score();

            cerr << j << ' ' << j_score << ' ' << j_score_unmerged << '\n';
            orientation_distribution.update(contact_graph);
        }

        path output_path = output_dir / ("orientations_" + to_string(i) + ".csv");
        orientation_distribution.write_contact_map(output_path, id_map);

        vector <pair <orientation_edge_t, orientation_weight_t> > ordered_edges;
        ordered_edges.reserve(contact_graph.edge_count());

        for (const auto& [edge,weights]: orientation_distribution.edge_weights){
            auto current_weight = max(weights[0],weights[1]);

            if (current_weight == sample_size){
                ordered_edges.emplace_back(edge, weights);
            }
        }

        sort(ordered_edges.begin(), ordered_edges.end(), [&](
                const pair<orientation_edge_t,orientation_weight_t>& a,
                const pair<orientation_edge_t,orientation_weight_t>& b
                ){

            auto a_ordinal = max(a.second[0], a.second[1]);
            auto b_ordinal = max(b.second[0], b.second[1]);

            bool result = a_ordinal > b_ordinal;

            if (a_ordinal == b_ordinal){
                auto a0_consistency = contact_graph.compute_consistency_score(a.first.first);
                auto a1_consistency = contact_graph.compute_consistency_score(a.first.second);
                auto b0_consistency = contact_graph.compute_consistency_score(b.first.first);
                auto b1_consistency = contact_graph.compute_consistency_score(b.first.second);
                auto a_avg = (a0_consistency + a1_consistency) / 2;
                auto b_avg = (b0_consistency + b1_consistency) / 2;

                result = a_avg > b_avg;
            }

            return result;
        });

        auto top_result = ordered_edges.front().second;
        auto max_weight = max(top_result[0],top_result[1]);
        auto current_weight = max_weight;

        alt_component_t component_a;
        alt_component_t component_b;

        unordered_set<int32_t> visited_nodes;

        // TODO: Test stability

//        cerr << "Starting with max_weight: " << max_weight << '\n';
        size_t e = 0;
        for (auto& [edge, weights]: ordered_edges){
            if (double(e)/double(ordered_edges.size()) > 0.2){
                break;
            }

            current_weight = max(weights[0],weights[1]);

            if (current_weight < sample_size){
                break;
            }

            auto partition = contact_graph.get_partition(edge.first);
            bool flipped = weights[0] < weights[1];

            if (visited_nodes.count(edge.first) + visited_nodes.count(edge.second) > 0){
                continue;
            }

            cerr << id_map.get_name(edge.first) << ' ' << id_map.get_name(edge.second) << ' ' << weights[0] << ' ' << weights[1] << ' ' << current_weight << '\n';

            contact_graph.get_alt_component(edge.first, false, component_a);
            contact_graph.get_alt_component(edge.second, false, component_b);

            if (flipped){
                flip_component(component_b);
            }

            cerr << "merging:" << '\n';
            cerr << "\ta" << '\n';
            for (auto& id: component_a.first){
                cerr << "\t0 " << id_map.get_name(id) << ' ' << int(contact_graph.get_partition(id)) << '\n';
            }
            for (auto& id: component_a.second){
                cerr << "\t1 " << id_map.get_name(id) << ' ' << int(contact_graph.get_partition(id)) << '\n';
            }
            cerr << "\tb" << '\n';
            for (auto& id: component_b.first){
                cerr << "\t0 " << id_map.get_name(id) << ' ' << int(contact_graph.get_partition(id)) << '\n';
            }
            for (auto& id: component_b.second){
                cerr << "\t1 " << id_map.get_name(id) << ' ' << int(contact_graph.get_partition(id)) << '\n';
            }

            contact_graph.add_alt(component_a, component_b, true);
            contact_graph.set_partition(edge.first, partition);

            alt_component_t component_merged;
            contact_graph.get_alt_component(edge.first, false, component_merged);
            cerr << "MERGED:" << '\n';
            for (auto& id: component_merged.first){
                cerr << "\t0 " << id_map.get_name(id) << ' ' << int(contact_graph.get_partition(id)) << '\n';
            }
            for (auto& id: component_merged.second){
                cerr << "\t1 " << id_map.get_name(id) << ' ' << int(contact_graph.get_partition(id)) << '\n';
            }

            for (auto& id: component_merged.first) {
                visited_nodes.emplace(id);
            }

            for (auto& id: component_merged.second) {
                visited_nodes.emplace(id);
            }
        }
    }

    cerr << "Final phase:" << '\n';
    auto best_score = numeric_limits<double>::min();
    auto best_phase_state = phase_state;

    for (size_t i=0; i<sample_size; i++){
        phase_contacts(contact_graph, n_threads, m_iterations*3);
        contact_graph.get_partitions(phase_state);
        unmerged_contact_graph.set_partitions(phase_state);

        auto score_unmerged = unmerged_contact_graph.compute_total_consistency_score();
        cerr << i << ' ' << score_unmerged << '\n';

        if (score_unmerged > best_score){
            best_score = score_unmerged;
            best_phase_state = phase_state;
        }
    }

    contact_graph.get_partitions(best_phase_state);

    cerr << "Final score:" << ' ' << best_score << '\n';
}


void phase_contacts(
        MultiContactGraph& contact_graph,
        size_t n_threads,
        size_t m_iterations,
        bool component_scale
){

    vector<thread> threads;
    vector <pair <int32_t,int8_t> > best_partitions;
    vector<int32_t> ids = {};
    atomic<double> best_score = std::numeric_limits<double>::min();
    atomic<size_t> job_index = 0;
    mutex phase_mutex;

    contact_graph.randomize_partitions();
    contact_graph.get_partitions(best_partitions);

//    cerr << "start score: " << contact_graph.compute_total_consistency_score() << '\n';

    auto f = random_node_scale_phase_search;

    if (component_scale){
        f = random_component_scale_phase_search;
        contact_graph.get_alt_component_representatives(ids);
    }
    else{
        contact_graph.get_node_ids(ids);
    }

    // Launch threads
    for (uint64_t i=0; i<n_threads; i++){
        try {
            threads.emplace_back(thread(
                    f,
                    contact_graph,
                    cref(ids),
                    ref(best_partitions),
                    ref(best_score),
                    ref(job_index),
                    ref(phase_mutex),
                    m_iterations
            ));
        } catch (const exception &e) {
            cerr << e.what() << "\n";
            exit(1);
        }
    }

    // Wait for threads to finish
    for (auto& t: threads){
        t.join();
    }

    contact_graph.set_partitions(best_partitions);
}


}