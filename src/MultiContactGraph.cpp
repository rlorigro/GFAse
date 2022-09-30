#include "MultiContactGraph.hpp"

#include <ostream>
#include <queue>

using std::ofstream;
using std::ostream;
using std::set_intersection;
using std::queue;
using std::cerr;
using std::min;
using std::max;


namespace gfase{


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


const array<string,3> MultiContactGraph::colors = {"Cornflower Blue", "Plum", "Tomato"};


MultiContactGraph::MultiContactGraph(const contact_map_t& contact_map, const IncrementalIdMap<string>& id_map){
    nodes.set_deleted_key(-1);

    string s;

    for (const auto& [a,sub_map]: contact_map){
        try_insert_node(a);

        for (const auto& [b,count]: sub_map){
            try_insert_node(b);
            insert_edge(a,b,count);
        }
    }
}


MultiContactGraph::MultiContactGraph(){
    nodes.set_deleted_key(-1);
}


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


void MultiContactGraph::add_alt(int32_t a, int32_t b){
    if (nodes.count(a) == 0 or nodes.count(b) == 0){
        throw runtime_error("ERROR: cannot add alt with nonexistent node id: (" + to_string(a) + "," + to_string(b) +")");
    }

    if (a == b){
        throw runtime_error("ERROR: cannot add alt to itself: " + to_string(b));
    }

    // Start by doing alt-wise BFS to get bipartite component of nodes
    alt_component_t component;
    get_alt_component(a, false, component);

    if (component.first.count(b) == 0){
        component.second.emplace(b);
    }
    else{
        for (auto& item: component.first){
            cerr << 0 << ' ' << item << '\n';
        }
        for (auto& item: component.second){
            cerr << 1 << ' ' << item << '\n';
        }

        throw runtime_error("ERROR: adding alt for " + to_string(a) + ',' + to_string(b) + " would result in non-bipartite component, see above for details");
    }
    if (component.second.count(a) == 0) {
        component.first.emplace(a);
    }
    else{
        for (auto& item: component.first){
            cerr << 0 << ' ' << item << '\n';
        }
        for (auto& item: component.second){
            cerr << 1 << ' ' << item << '\n';
        }

        throw runtime_error("ERROR: adding alt for " + to_string(a) + ',' + to_string(b) + " would result in non-bipartite component, see above for details");
    }

    auto& node_b = nodes.at(b);
    auto& node_a = nodes.at(a);

    for (auto& id: component.first){
        auto& node_n = nodes.at(id);

        // Enforce all-vs-all connectivity in alt components
        node_n.alts.emplace(b);
        node_b.alts.emplace(id);

        // No valid edge can exist between alts
        remove_edge(id, b);

        // Assign arbitrary partition
        node_n.partition = -1;
    }
    for (auto& id: component.second){
        auto& node_n = nodes.at(id);

        // Enforce all-vs-all connectivity in alt components
        node_n.alts.emplace(a);
        node_a.alts.emplace(id);

        // No valid edge can exist between alts
        remove_edge(a, id);

        // Assign arbitrary partition
        node_n.partition = 1;
    }
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


size_t MultiContactGraph::size(){
    return nodes.size();
}


void MultiContactGraph::resize(){
    nodes.resize(0);
}


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


double MultiContactGraph::compute_total_consistency_score() const{
    double score = 0;

    for_each_edge([&](const pair<int32_t,int32_t> edge, int32_t weight){
        auto& a = nodes.at(edge.first);
        auto& b = nodes.at(edge.second);
        score += get_score(a, b, weight);
    });

    return score;
}


void MultiContactGraph::get_partitions(vector <pair <int32_t,int8_t> >& partitions) const{
    partitions.clear();
    partitions.resize(nodes.size());

    size_t i = 0;
    for (auto& [n, node]: nodes){
        partitions[i] = {n,node.partition};
        i++;
    }
}


void MultiContactGraph::randomize_partitions(){
    // True random number
    std::random_device rd;

    // Pseudorandom generator with true random seed
    std::mt19937 rng(rd());
    std::uniform_int_distribution<int> uniform_distribution(0,2);

    for (auto& [id,node]: nodes){
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


void MultiContactGraph::set_partitions(const vector <pair <int32_t,int8_t> >& partitions){
    for (auto& [n, p]: partitions){
        set_partition(n, p);
    }
}


void MultiContactGraph::get_node_ids(vector<int32_t>& ids){
    ids.resize(nodes.size());

    size_t i = 0;

    for (auto& [id,node]: nodes){
        ids[i] = id;
        i++;
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


void MultiContactGraph::write_bandage_csv(path output_path, IncrementalIdMap<string>& id_map) const{
    ofstream file(output_path);

    if (not file.is_open() or not file.good()) {
        throw std::runtime_error("ERROR: could not write to file: " + output_path.string());
    }

    file << "Name" << ',' << "Phase" << ',' << "Coverage" << ',' << "Length" << ',' << "Color" << '\n';

    for (auto& [id,node]: nodes){
        auto name = id_map.get_name(id);
        file << name << ',' << int(node.partition) << ',' << node.coverage << ',' << node.length << ',' << colors[node.partition+1] << '\n';
    }
}


void MultiContactGraph::write_node_data(path output_path, IncrementalIdMap<string>& id_map) const{
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


void random_multicontact_phase_search(
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
    std::uniform_int_distribution<int> uniform_distribution(0,int(contact_graph.size()-1));

    contact_graph.set_partitions(best_partitions);

    double total_score;

    while (m < m_iterations) {
        // Randomly perturb
        for (size_t i=0; i<((contact_graph.size()/30) + 1); i++) {
            auto r = ids[uniform_distribution(rng)];

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

        for (size_t i=0; i<contact_graph.size(); i++) {
            auto n = ids[uniform_distribution(rng)];

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

        cerr << m << ' ' << best_score << ' ' << total_score << ' ' << std::flush;
//        for (auto& [n,p]: best_partitions){
//            cerr << '(' << n << ',' << int(p) << ") ";
//        }
        cerr << '\n';

        phase_mutex.unlock();

        m = job_index.fetch_add(1);
    }

    contact_graph.set_partitions(best_partitions);
}



}