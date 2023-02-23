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
using std::ref;

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


void OrientationDistribution::update(const VectorMultiContactGraph& contact_graph){
    contact_graph.for_each_edge([&](const pair<int32_t,int32_t> edge, int32_t weight){
        auto [a,b] = edge;

        bool orientation = contact_graph.get_partition(a) == contact_graph.get_partition(b);

        edge_weights[edge][orientation]++;
    });
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


void random_phase_search(VectorMultiContactGraph& contact_graph, size_t m_iterations){
    vector <pair <int32_t,int8_t> > best_partitions;
    double best_score = std::numeric_limits<double>::min();

    vector<int32_t> ids = {};
    contact_graph.get_node_ids(ids);

    contact_graph.randomize_partitions();
    contact_graph.get_partitions(best_partitions);

    // True random number
    std::random_device rd;

    // Pseudorandom generator with true random seed
    std::mt19937 rng(rd());
    std::uniform_int_distribution<int> uniform_distribution(0,int(ids.size()-1));

    double total_score;

    for (size_t m=0; m<m_iterations; m++) {
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

        for (size_t i=0; i<ids.size()*3; i++) {
            auto n = ids.at(uniform_distribution(rng));

            if (contact_graph.edge_count(n) == 0){
                continue;
            }
            bool has_alt = contact_graph.has_alt(n);
            auto prev_score = contact_graph.compute_consistency_score(n);
            auto prev_partition = contact_graph.get_partition(n);

            double max_score = prev_score;
            int8_t p_max = prev_partition;

            if (prev_partition == -1 or prev_partition == 0) {
//                cerr << "TEST P = 1" << '\n';

                auto score = contact_graph.compute_consistency_score(n, 1);

                if (score > max_score) {
                    max_score = score;
                    p_max = 1;
                }
            }

            if (prev_partition == 1 or prev_partition == 0) {
//                cerr << "TEST P = -1" << '\n';

                auto score = contact_graph.compute_consistency_score(n, -1);

                if (score > max_score) {
                    max_score = score;
                    p_max = -1;
                }
            }

            // If the node has no "alt" it can be made neutral
            if ((not has_alt) and prev_partition != 0){
                auto score = contact_graph.compute_consistency_score(n, 0);

                if (score > max_score) {
                    max_score = score;
                    p_max = 0;
                }
            }

            contact_graph.set_partition(n, p_max);
        }

        total_score = contact_graph.compute_total_consistency_score();

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

    }
}


void flip_component(alt_component_t& c){
    auto temp = c.second;
    c.second = c.first;
    c.first = temp;
}


void sample_with_threads(vector<VectorMultiContactGraph>& contact_graphs_per_thread,
                         size_t core_iterations,
                         atomic<size_t>& job_index){
    auto i = job_index.fetch_add(1);

    while (i < contact_graphs_per_thread.size()){
        random_phase_search(contact_graphs_per_thread[i], core_iterations);
        i = job_index.fetch_add(1);
    }
}


void sample_orientation_distribution(
        OrientationDistribution& orientation_distribution,
        MultiContactGraph& contact_graph,
        size_t sample_size,
        size_t n_threads,
        size_t core_iterations
        ){

    vector<thread> threads;
    vector<int32_t> ids = {};

    contact_graph.randomize_partitions();
    contact_graph.get_node_ids(ids);

    vector<VectorMultiContactGraph> contact_graphs_per_thread(sample_size,contact_graph);
    atomic<size_t> job_index = 0;

    // Launch threads
    for (uint64_t i=0; i<n_threads; i++){
        try {
            threads.emplace_back(thread(
                    sample_with_threads,
                    ref(contact_graphs_per_thread),
                    core_iterations,
                    ref(job_index)
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

    double best_score = numeric_limits<double>::min();
    vector <pair <int32_t,int8_t> > best_partitions;

    cerr << "sampling results: " << '\n';
    for (const auto& result: contact_graphs_per_thread){
        auto score = result.compute_total_consistency_score();

        if (score > best_score){
            best_score = score;
            result.get_partitions(best_partitions);
        }

        cerr << score << '\n';
        orientation_distribution.update(result);
    }

    contact_graph.set_partitions(best_partitions);
}


void monte_carlo_phase_contacts(
        MultiContactGraph& contact_graph,
        const IncrementalIdMap<string>& id_map,
        size_t core_iterations,
        size_t sample_size,
        size_t n_rounds,
        size_t n_threads,
        path output_dir
        ){

    // Keep the original graph for scoring purposes (some bubbles will be merged later)
    MultiContactGraph unmerged_contact_graph = contact_graph;

    // Don't evaluate a pair of bubbles twice during merging
    unordered_set<orientation_edge_t> visited_edges;
    vector <pair <int32_t,int8_t> > phase_state;

    for (size_t i=0; i<n_rounds; i++){
        // Initialize DS for tracking results of repeated samples from the converged graph
        OrientationDistribution orientation_distribution(contact_graph);

        cerr << "---- " << i << " ----" << '\n';
        sample_orientation_distribution(
                orientation_distribution,
                contact_graph,
                sample_size,
                n_threads,
                core_iterations);

        // Convert to non-mutable graph for efficiency of optimization
        VectorMultiContactGraph vector_contact_graph(contact_graph);

        path components_path = output_dir / ("components_" + to_string(i) + ".csv");
        vector_contact_graph.write_alt_components(components_path, id_map);

        path orientations_path = output_dir / ("orientations_" + to_string(i) + ".csv");
        orientation_distribution.write_contact_map(orientations_path, id_map);

        // Initialize storage for the edges in order of best first (for merging purposes)
        vector <pair <orientation_edge_t, orientation_weight_t> > ordered_edges;
        ordered_edges.reserve(contact_graph.edge_count());

        // Only accumulate edges which are perfectly consistent
        for (const auto& [edge,weights]: orientation_distribution.edge_weights){
            auto current_weight = max(weights[0],weights[1]);

            if (current_weight == sample_size){
                ordered_edges.emplace_back(edge, weights);
            }
        }

        // Order by edge weight
        sort(ordered_edges.begin(), ordered_edges.end(), [&](
                const pair<orientation_edge_t,orientation_weight_t>& a,
                const pair<orientation_edge_t,orientation_weight_t>& b
                ){

            auto a_ordinal = max(a.second[0], a.second[1]);
            auto b_ordinal = max(b.second[0], b.second[1]);

            bool result = a_ordinal > b_ordinal;

            if (a_ordinal == b_ordinal){
                auto a0_consistency = vector_contact_graph.compute_consistency_score(a.first.first);
                auto a1_consistency = vector_contact_graph.compute_consistency_score(a.first.second);
                auto b0_consistency = vector_contact_graph.compute_consistency_score(b.first.first);
                auto b1_consistency = vector_contact_graph.compute_consistency_score(b.first.second);
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

        // Iterate top 20% of edges
        size_t e = 0;
        for (auto& [edge, weights]: ordered_edges){
            if (double(e)/double(ordered_edges.size()) > 0.2){
                break;
            }

            current_weight = max(weights[0],weights[1]);

            if (current_weight < sample_size){
                break;
            }

            // Keep track of the orientation so that merging step merges in correct orientation
            auto partition = vector_contact_graph.get_partition(edge.first);
            bool flipped = weights[0] < weights[1];

            if (visited_nodes.count(edge.first) + visited_nodes.count(edge.second) > 0){
                continue;
            }

            vector_contact_graph.get_alt_component(edge.first, false, component_a);
            vector_contact_graph.get_alt_component(edge.second, false, component_b);

            if (flipped){
                flip_component(component_b);
            }

            // Update alt relationships among nodes within merged bubbles, and set partitions
            contact_graph.add_alt(component_a, component_b, true);
            contact_graph.set_partition(edge.first, partition);

            alt_component_t component_merged;
            contact_graph.get_alt_component(edge.first, false, component_merged);

            for (auto& id: component_merged.first) {
                visited_nodes.emplace(id);
            }

            for (auto& id: component_merged.second) {
                visited_nodes.emplace(id);
            }
        }
    }

    OrientationDistribution orientation_distribution(contact_graph);

    // Perform finishing convergence on the most merged graph, with more iterations
    cerr << "Final phase:" << '\n';
    sample_orientation_distribution(
            orientation_distribution,
            contact_graph,
            sample_size,
            n_threads,
            3*core_iterations);

    // Store best result for future use
    vector <pair <int32_t,int8_t> > best_partitions;
    contact_graph.get_partitions(best_partitions);
    unmerged_contact_graph.set_partitions(best_partitions);

    auto final_unmerged_score = unmerged_contact_graph.compute_total_consistency_score();
    cerr << "Final unmerged score: " << final_unmerged_score << '\n';

    VectorMultiContactGraph g(contact_graph);

    path components_path = output_dir / ("components_final.csv");
    g.write_alt_components(components_path, id_map);

    path orientations_path = output_dir / ("orientations_final.csv");
    orientation_distribution.write_contact_map(orientations_path, id_map);

    // Reset the contact graph to the unmerged state so its alts can be used for chaining in future methods
    contact_graph = unmerged_contact_graph;
}


}
