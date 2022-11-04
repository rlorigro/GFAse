#include "VectorMultiContactGraph.hpp"

#include <queue>

using std::queue;


namespace gfase{


VectorMultiNode::VectorMultiNode(int8_t partition):
        neighbors(),
        coverage(0),
        length(0),
        partition(partition)
{}


VectorMultiNode::VectorMultiNode(const MultiNode& node):
        neighbors(node.neighbors.size()),
        alts(node.alts.size()),
        coverage(node.coverage),
        length(node.length),
        partition(node.partition)
{
    size_t i;

    i = 0;
    for (auto id: node.neighbors){
        neighbors[i] = {id,0};
        i++;
    }

    i = 0;
    for (auto id: node.alts){
        alts[i] = id;
        i++;
    }
}


bool VectorMultiNode::has_alt() const{
    return not alts.empty();
}


VectorMultiContactGraph::VectorMultiContactGraph(const MultiContactGraph& contact_graph):
        nodes(contact_graph.get_max_id()+1, nullptr),    // Fill with -1 for anticipated empty values (gaps)
        edge_weights(contact_graph.get_max_id()+1)
{
    contact_graph.for_each_node([&](int32_t id, const MultiNode& n){
        VectorMultiNode vn(n);
        nodes.at(id) = make_shared<VectorMultiNode>(vn);
    });

    contact_graph.for_each_edge([&](const pair<int32_t,int32_t> edge, int32_t weight){
        edge_weights[edge.first].emplace(edge.second, weight);
        for (auto& neighbor: nodes.at(edge.first)->neighbors){
            if (neighbor.first == edge.second){
                neighbor.second = weight;
            }
        }
    });
}


/// Use BFS on node alts to get connected component that represents a bubble
/// \param id
/// \param validate
/// \param component
void VectorMultiContactGraph::get_alt_component(int32_t id, bool validate, alt_component_t& component) const{
    component = {};

    queue <pair <int32_t,int32_t> > q;
    q.emplace(id,0);

    while(not q.empty()){
        auto& [current_id,distance] = q.front();
        q.pop();

        auto& node = nodes.at(current_id);

        if (distance % 2 == 0){
            component.first.emplace(current_id);
        }
        else{
            component.second.emplace(current_id);
        }

        for (auto& alt_id: node->alts) {
            auto a_result = component.first.find(alt_id);
            auto b_result = component.second.find(alt_id);

            bool a_found = (a_result != component.first.end());
            bool b_found = (b_result != component.second.end());

            if (not (a_found or b_found)){
                q.emplace(alt_id, distance+1);
            }
        }
    }
}


void VectorMultiContactGraph::set_partition(int32_t id, int8_t partition) {
    auto& node = nodes.at(id);

    node->partition = partition;

    // If this node is linked to an alt, alt must be maintained in an opposite state,
    // and double-alts must be maintained in identical state
    if (node->has_alt()){
        if (partition == 0) {
            throw runtime_error("ERROR: cannot set 0 partition for bubble: " + to_string(id));
        }

        alt_component_t component;
        get_alt_component(id, false, component);

        for (auto& alt_id: component.second){
            auto& alt = nodes.at(alt_id);
            alt->partition = partition*int8_t(-1);
        }
        for (auto& alt_id: component.first){
            auto& alt = nodes.at(alt_id);
            alt->partition = partition;
        }
    }
}


double VectorMultiContactGraph::get_score(const shared_ptr<VectorMultiNode>& a, const shared_ptr<VectorMultiNode>& b, int32_t weight){
    double score = 0;

    auto p_a = a->partition;
    auto p_b = b->partition;

    if (p_a != 0 and p_b != 0) {
        score = p_a * p_b * weight;
    }
    else{
        score = 0;
    }

    return score;
}


double VectorMultiContactGraph::compute_consistency_score(int32_t id) const{
    double score = 0;

    const auto n = nodes.at(id);

    if (n == nullptr){
        throw runtime_error("ERROR: VectorMultiContactGraph::compute_consistency_score: nonexistent node ID: " + to_string(id));
    }

    for (auto& [id_other, weight]: n->neighbors) {
        // Skip self edges if there are any
        if (id == id_other) {
            continue;
        }

        const auto& n_other = nodes.at(id_other);

        score += get_score(n, n_other, weight);
//        cerr << '\t' << id << "<->" << id_other << ' ' << int(n.partition) << 'x' << int(n_other.partition) << 'x' << edge_weights.at(edge(id, id_other)) << ' ' << score << '\n';
    }

    for (auto alt_id: n->alts){
        const auto& n_alt = nodes.at(alt_id);

        for (auto& [id_other,weight]: n_alt->neighbors) {
            auto& n_other = nodes.at(id_other);

            score += get_score(n_alt, n_other, weight);
        }
    }

    return score;
}


double VectorMultiContactGraph::compute_total_consistency_score() const{
    double score = 0;

    for (int32_t id_a=0; id_a<edge_weights.size(); id_a++){
        if (nodes.at(id_a) == nullptr){
            continue;
        }

        for (auto& [id_b, weight]: edge_weights[id_a]) {
            // Skip self edges if any exist
            if (id_a == id_b) {
                continue;
            }

            auto& a = nodes.at(id_a);
            auto& b = nodes.at(id_b);
            score += get_score(a, b, weight);
        }
    }

    return score;
}


void VectorMultiContactGraph::get_partitions(vector <pair <int32_t,int8_t> >& partitions) const{
    partitions.clear();

    for (size_t n=0; n<nodes.size(); n++){
        auto& node = nodes.at(n);
        if (node != nullptr) {
            partitions.emplace_back(n,node->partition);
        }
    }
}


void VectorMultiContactGraph::set_partitions(const vector <pair <int32_t,int8_t> >& partitions){
    for (const auto& [n, p]: partitions){
        set_partition(n, p);
    }
}


size_t VectorMultiContactGraph::edge_count(int32_t id) const{
    return nodes.at(id)->neighbors.size();
}


bool VectorMultiContactGraph::has_alt(int32_t id) const{
    return nodes.at(id)->has_alt();
}


int8_t VectorMultiContactGraph::get_partition(int32_t id) const{
    return nodes.at(id)->partition;
}


}