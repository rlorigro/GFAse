#ifndef GFASE_PHASE_TRIPARTITION_HPP
#define GFASE_PHASE_TRIPARTITION_HPP

#include "VectorMultiContactGraph.hpp"
#include "MultiContactGraph.hpp"


namespace gfase{


using orientation_edge_t = pair <int32_t,int32_t>;
using orientation_weight_t = array<int32_t, 2>;


class OrientationDistribution{
public:
    unordered_map <orientation_edge_t, orientation_weight_t> edge_weights;
    vector <alt_component_t> alt_components;

    OrientationDistribution(const MultiContactGraph& contact_graph);
    void write_contact_map(path output_path, const IncrementalIdMap<string>& id_map) const;
    void update(const VectorMultiContactGraph& contact_graph);
    void update(const MultiContactGraph& contact_graph);
};


void random_phase_search(VectorMultiContactGraph& contact_graph, size_t m_iterations);


void sample_with_threads(
        vector<VectorMultiContactGraph>& contact_graphs_per_thread,
        size_t m_iterations,
        atomic<size_t>& job_index);


void sample_orientation_distribution(
        OrientationDistribution& orientationDistribution,
        MultiContactGraph& contact_graph,
        size_t sample_size,
        size_t n_threads,
        size_t m_iterations
);


void monte_carlo_phase_contacts(
        MultiContactGraph& contact_graph,
        const IncrementalIdMap<string>& id_map,
        path output_dir,
        size_t n_threads
);


}

#endif //GFASE_PHASE_TRIPARTITION_HPP
