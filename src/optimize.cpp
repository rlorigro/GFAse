#include "optimize.hpp"
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



//void random_component_scale_phase_search(
//        MultiContactGraph& contact_graph,
//        const vector<int32_t>& ids,
//        vector <pair <int32_t,int8_t> >& best_partitions,
//        atomic<double>& best_score,
//        atomic<size_t>& job_index,
//        mutex& phase_mutex,
//        size_t m_iterations
//){
//
//    size_t m = job_index.fetch_add(1);
//
//    // True random number
//    std::random_device rd;
//
//    // Pseudorandom generator with true random seed
//    std::mt19937 rng(rd());
//    std::uniform_int_distribution<int> uniform_distribution(0,int(ids.size()-1));
//
//    contact_graph.set_partitions(best_partitions);
//
//    alt_component_t component;
//
//    double total_score;
//
//    while (m < m_iterations) {
//        // Randomly perturb
//        for (size_t i=0; i<((ids.size()/30) + 1); i++) {
//            auto r = ids.at(uniform_distribution(rng));
//
//            int8_t p;
//            if (contact_graph.has_alt(r)){
//                // Only allow {1,-1}
//                p = int8_t((uniform_distribution(rng) % 2));
//
//                if (p == 0){
//                    p = -1;
//                }
//            }
//            else{
//                // Allow {1,0,-1}
//                p = int8_t((uniform_distribution(rng) % 3) - 1);
//            }
//
//            contact_graph.set_partition(r, p);
//        }
//
//        for (size_t i=0; i<ids.size(); i++) {
//            auto n = ids.at(uniform_distribution(rng));
//
//            int64_t max_score = std::numeric_limits<int64_t>::min();
//            int8_t p_max = 0;
//
//            contact_graph.get_alt_component(n, false, component);
//            bool has_alt = not component.first.empty() and not component.second.empty();
//
//            for (int8_t p=-1; p<=1; p++) {
//                // If the node has an "alt" it can't be made neutral
//                if (has_alt and p==0){
//                    continue;
//                }
//
//                contact_graph.set_partition(component, p);
//
//                auto score = contact_graph.compute_consistency_score(component);
//
//                if (score > max_score) {
//                    max_score = score;
//                    p_max = p;
//                }
//
////                cerr << n << ' ' << int(has_alt) << ' ' << int(p) << ' ' << score << ' ' << max_score << '\n';
//            }
//
//            contact_graph.set_partition(n, p_max);
//        }
//
//        total_score = contact_graph.compute_total_consistency_score();
//
//        phase_mutex.lock();
//        if (total_score > best_score) {
//            best_score = total_score;
//            contact_graph.get_partitions(best_partitions);
//        }
//        else {
//            contact_graph.set_partitions(best_partitions);
//        }
//
////        cerr << m << ' ' << best_score << ' ' << total_score << ' ' << std::flush;
////        for (auto& [n,p]: best_partitions){
////            cerr << '(' << n << ',' << int(p) << ") ";
////        }
////        cerr << '\n';
//
//        phase_mutex.unlock();
//
//        m = job_index.fetch_add(1);
//    }
//
//    contact_graph.set_partitions(best_partitions);
//}


void random_node_scale_phase_search(
        VectorMultiContactGraph contact_graph,
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

    while (job_index < m_iterations) {
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

            if (contact_graph.edge_count(n) == 0){
                continue;
            }
//            cerr << "---- " << n << " ----" << '\n';
            bool has_alt = contact_graph.has_alt(n);
            auto prev_score = contact_graph.compute_consistency_score(n);
            auto prev_partition = contact_graph.get_partition(n);

            double max_score = prev_score;
            int8_t p_max = prev_partition;

            // TODO: don't actually "set" partition, pass values to be used
            // TODO: find convergence issue
            // TODO: cache components?

//            cerr << "max_score: " << ' ' << max_score << '\n';
//            cerr << "p_max: " << ' ' << int(p_max) << '\n';
//            cerr << "prev_score: " << ' ' << prev_score << '\n';
//            cerr << "prev_partition: " << ' ' << int(prev_partition) << '\n';

            if (prev_partition == -1 or prev_partition == 0) {
//                cerr << "TEST P = 1" << '\n';

                contact_graph.set_partition(n, 1);
                auto score = contact_graph.compute_consistency_score(n);

                if (score > max_score) {
                    max_score = score;
                    p_max = 1;
                }

//                cerr << "max_score: " << ' ' << max_score << '\n';
//                cerr << "score: " << ' ' << score << '\n';
//                cerr << "p_max: " << ' ' << int(p_max) << '\n';

            }

            if (prev_partition == 1 or prev_partition == 0) {
//                cerr << "TEST P = -1" << '\n';

                contact_graph.set_partition(n, -1);
                auto score = contact_graph.compute_consistency_score(n);

                if (score > max_score) {
                    max_score = score;
                    p_max = -1;
                }

//                cerr << "max_score: " << ' ' << max_score << '\n';
//                cerr << "score: " << ' ' << score << '\n';
//                cerr << "p_max: " << ' ' << int(p_max) << '\n';

            }

            // If the node has no "alt" it can be made neutral
            if ((not has_alt) and prev_partition != 0){
                contact_graph.set_partition(n, 0);
                auto score = contact_graph.compute_consistency_score(n);

                if (score > max_score) {
                    max_score = score;
                    p_max = 0;
                }
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

//        cerr << m << ' ' << best_score << ' ' << total_score << ' ';
//        size_t x = 0;
//        for (auto& [n,p]: best_partitions){
//            if (++x > 10){
//                break;
//            }
//            cerr << '(' << n << ',' << int(p) << ") ";
//        }
//        cerr << '\n';

        phase_mutex.unlock();

        m = job_index.fetch_add(1);
    }
}


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


void flip_component(alt_component_t& c){
    auto temp = c.second;
    c.second = c.first;
    c.first = temp;
}


void phase_contacts(
        MultiContactGraph& contact_graph,
        size_t n_threads,
        size_t m_iterations
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

//    if (component_scale){
//        f = random_component_scale_phase_search;
//        contact_graph.get_alt_component_representatives(ids);
//    }
//    else{
        contact_graph.get_node_ids(ids);
//    }

    VectorMultiContactGraph vector_contact_graph(contact_graph);

    // Launch threads
    for (uint64_t i=0; i<n_threads; i++){
        try {
            threads.emplace_back(thread(
                    f,
                    vector_contact_graph,
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


}