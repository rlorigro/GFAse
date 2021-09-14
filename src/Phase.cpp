#include "FixedBinarySequence.hpp"
#include "Phase.hpp"

namespace gfase{


// Find any nodes that are adjacent to the beginning and end of a path, as long as they are the only adjacent node
pair<handle_t, bool> find_singleton_adjacent_handle(const PathHandleGraph& graph, const handle_t& h, bool left) {
    handle_t adjacent_handle;
    size_t n_adjacent = 0;
    bool success = false;

    graph.follow_edges(h, left, [&](const handle_t& other_handle) {
        if (n_adjacent > 0) {
            return false;
        }

        adjacent_handle = other_handle;
        n_adjacent++;

        return true;
    });

    // Only return true if there was exactly one adjacent handle
    if (n_adjacent == 1) {
        success = true;
    }

    return {adjacent_handle, success};
}


void extend_paths(MutablePathMutableHandleGraph& graph) {
    vector<pair<path_handle_t, handle_t> > to_be_prepended;
    vector<pair<path_handle_t, handle_t> > to_be_appended;

    graph.for_each_path_handle([&](const path_handle_t& p) {
        auto begin_handle = graph.get_handle_of_step(graph.path_begin(p));
        auto end_handle = graph.get_handle_of_step(graph.path_back(p));

        auto left_result = find_singleton_adjacent_handle(graph, begin_handle, true);
        auto right_result = find_singleton_adjacent_handle(graph, end_handle, false);

        if (left_result.second) {
            size_t n_paths = 0;
            graph.for_each_step_on_handle(left_result.first, [&](const step_handle_t& s) {
                n_paths++;
            });

            // Verify that there is not already some other path covering the adjacent node
            if (n_paths == 0) {
                to_be_prepended.emplace_back(p, left_result.first);
            } else {
                throw runtime_error(
                        "ERROR: node left of haplotype path " + graph.get_path_name(p) + " is covered by other paths");
            }
        }

        if (right_result.second) {
            size_t n_paths = 0;
            graph.for_each_step_on_handle(right_result.first, [&](const step_handle_t& s) {
                n_paths++;
            });

            // Verify that there is not already some other path covering the adjacent node
            if (n_paths == 0) {
                to_be_appended.emplace_back(p, right_result.first);
            } else {
                throw runtime_error(
                        "ERROR: node right haplotype path " + graph.get_path_name(p) + " is covered by other paths");
            }
        }
    });

    for (auto& item: to_be_prepended) {
        // PREpend the LEFT side node if it meets the conditions
        graph.prepend_step(item.first, item.second);
    }

    for (auto& item: to_be_appended) {
        // Append the RIGHT side node if it meets the conditions
        graph.append_step(item.first, item.second);
    }
}


void phase_haplotype_paths(path gfa_path, size_t k, path paternal_kmers, path maternal_kmers) {
    HashGraph graph;
    IncrementalIdMap<string> id_map;
    KmerSets<FixedBinarySequence<uint64_t, 2> > ks(paternal_kmers, maternal_kmers);

    gfa_to_handle_graph(graph, id_map, gfa_path);

    // plot_graph(graph, "start_graph");

    cerr << "Extending paths into haploid regions..." << '\n';

    extend_paths(graph);

    cerr << "Iterating path kmers..." << '\n';

    cerr << "Number of components in graph: " << graph.get_path_count() << '\n';

    bool prev_has_diploid;

    // Iterate paths and for each node, collect kmers if node is only covered by one path
    graph.for_each_path_handle([&](const path_handle_t& p) {
        cerr << ">" << graph.get_path_name(p) << '\n';

        HaplotypePathKmer kmer(graph, p, k);

        string component_path_string = graph.get_path_name(p);

        while (kmer.step()) {
            if (kmer.has_diploid) {
                FixedBinarySequence<uint64_t,2> s(kmer.sequence);

                // compare kmer to parental kmers
                ks.increment_parental_kmer_count(component_path_string, s);

            } else {
                auto h = graph.get_handle_of_step(kmer.steps.back());
                auto length = graph.get_length(h);

                if (length > k + 1) {
                    kmer.initialize(kmer.steps.back(), length - k + 1);
                    prev_has_diploid = false;
                }
            }

            prev_has_diploid = kmer.has_diploid;
        }
    });
//    ks.normalize_kmer_counts();
    ks.print_component_parent_conf_matrix();
}

}