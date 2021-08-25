#include "IncrementalIdMap.hpp"
#include "gfa_to_handle.hpp"
#include "handle_to_gfa.hpp"
#include "GraphUtility.hpp"
#include "Filesystem.hpp"
#include "GfaReader.hpp"
#include "CLI11.hpp"

#include "bdsg/hash_graph.hpp"

#include <string>

using gfase::IncrementalIdMap;
using gfase::for_each_connected_component;
using gfase::split_connected_components;
using gfase::handle_graph_to_gfa;
using gfase::print_graph_paths;
using gfase::plot_graph;
using ghc::filesystem::path;

using bdsg::HashGraph;
using bdsg::MutablePathMutableHandleGraph;
using bdsg::MutablePathDeletableHandleGraph;
using handlegraph::path_handle_t;
using handlegraph::step_handle_t;
using handlegraph::handle_t;

using std::string;
using std::cout;
using std::cerr;


bool is_haplotype_bubble(const PathHandleGraph& graph, step_handle_t s){
    auto h = graph.get_handle_of_step(s);

    size_t n_other_paths_on_handle = 0;
    graph.for_each_step_on_handle(h, [&](const step_handle_t& s_other){
        if (s_other != s){
            n_other_paths_on_handle++;
        }
    });

    bool result;
    if (n_other_paths_on_handle == 0){
        // Must be haplotype bubble if paths are correct
        result = true;
    }
    else if (n_other_paths_on_handle == 1){
        // Not part of a bubble
        result = false;
    }
    else{
        throw runtime_error("ERROR: multiploid step in path for node: " + to_string(graph.get_id(h)));
        return result;
    }

    return result;
}


class HaplotypePathKmer{
public:
    /// Attributes ///
    const PathHandleGraph& graph;
    deque<step_handle_t> steps;
    deque<size_t> lengths;
    deque<bool> is_diploid;
    bool has_diploid;
    path_handle_t path;

    size_t start_index;
    size_t stop_index;

    step_handle_t terminal_step;
    size_t terminal_index;

    // K can be any size, so bit shift operations are not a simple option anymore
    deque<char> sequence;
    size_t k;

    /// Methods ///
    HaplotypePathKmer(const PathHandleGraph& graph,const path_handle_t& path, size_t k);

    // skip - jump to a position (step_handle_t, size_t) in the path and update internal records
    void initialize(step_handle_t s, size_t index);

    // step - walk a single bp forward and update internal records
    bool step();

    bool update_has_diploid();
};


HaplotypePathKmer::HaplotypePathKmer(const PathHandleGraph& graph, const path_handle_t& path, size_t k):
    graph(graph),
    has_diploid(false),
    path(path),
    k(k)
{
    if (k < 1){
        throw runtime_error("ERROR: k must be at least 1");
    }

    auto first_step = graph.path_begin(path);
    initialize(first_step, 0);
}


void HaplotypePathKmer::initialize(step_handle_t s, size_t index){
    auto h = graph.get_handle_of_step(s);
    bool step_is_diploid = is_haplotype_bubble(graph, s);

    steps.clear();
    lengths.clear();
    is_diploid.clear();
    sequence.clear();

    steps.emplace_back(s);
    lengths.emplace_back(graph.get_length(h));
    is_diploid.emplace_back(step_is_diploid);
    sequence.emplace_back(graph.get_base(h, index));
    update_has_diploid();

    start_index = index;
    stop_index = index;

    terminal_step = graph.path_back(path);
    terminal_index = graph.get_length(graph.get_handle_of_step(terminal_step));

    bool sufficient_path_length = true;

    while (sequence.size() < k - 1){
        // Preload the kmer queue up until it almost has k elements
        this->step();

        // Verify that this kmer actually fits in the path
        if (steps.back() == terminal_step and stop_index + 1 >= terminal_index) {
            sufficient_path_length = false;
        }
    }

    if (not sufficient_path_length){
        throw runtime_error("ERROR: path " + graph.get_path_name(path) + " does not have sufficient length for kmer size: " + to_string(k));
    }

}


// TODO: switch to standard queue (only need to pop front and push back)
bool HaplotypePathKmer::step(){
    bool has_next_step = true;

    // If it is safe to increment this node
    if (stop_index + 1 < lengths.back()){
        stop_index++;

        auto h = graph.get_handle_of_step(steps.back());
        sequence.emplace_back(graph.get_base(h,stop_index));

        // Use the deque like a cyclic queue
        if (sequence.size() > k) {
            sequence.pop_front();
        }
    }
    // If the iterator is at the end of the node
    else {
        if (steps.back() != terminal_step){
            step_handle_t next_step = graph.get_next_step(steps.back());
            auto next_handle = graph.get_handle_of_step(next_step);

            // Skip steps in path that are empty nodes
            while(graph.get_length(next_handle) == 0){
                next_step = graph.get_next_step(next_step);
                next_handle = graph.get_handle_of_step(next_step);
            }

            bool step_is_diploid = is_haplotype_bubble(graph, next_step);

            steps.emplace_back(next_step);
            lengths.emplace_back(graph.get_length(next_handle));
            is_diploid.emplace_back(step_is_diploid);
            sequence.emplace_back(graph.get_base(next_handle, 0));
            update_has_diploid();

            // Use the deque like a cyclic queue
            if (sequence.size() > k) {
                sequence.pop_front();
            }

            stop_index = 0;
        }
        else if (stop_index == lengths.back()){
            stop_index++;
        }
        else{
            return false;
        }
    }

    // Handle the trailing end of the kmer (queue front)
    if (start_index < lengths.front()) {
        if (sequence.size() == k) {
            start_index++;
        }
    }
    else{
        start_index = 0;
        lengths.pop_front();
        steps.pop_front();
        is_diploid.pop_front();
        update_has_diploid();
    }

    if (steps.back() == terminal_step and stop_index >= terminal_index){
        has_next_step = false;
    }

    return has_next_step;
}


bool HaplotypePathKmer::update_has_diploid(){
    bool has_diploid_nodes = false;
    for (size_t i=0; i<this->steps.size(); i++) {
        if (this->is_diploid[i]){
            has_diploid_nodes = true;
        }
    }

    has_diploid = has_diploid_nodes;
    return has_diploid_nodes;
}


void extract_haplotype_kmers_from_gfa(path gfa_path, size_t k){
    HashGraph graph;
    IncrementalIdMap<string> id_map;

    gfa_to_handle_graph(graph, id_map, gfa_path);

    cerr << "Iterating paths" << '\n';

    plot_graph(graph, "start_graph");

    bool prev_has_diploid;

    // Iterate paths and for each node, collect kmers if node is only covered by one path
    graph.for_each_path_handle([&](const path_handle_t& p){
        cerr << ">" << graph.get_path_name(p) << '\n';

        HaplotypePathKmer kmer(graph, p, k);

        while (kmer.step()){
            if (kmer.has_diploid){
                for (auto& c: kmer.sequence){
                    cerr << c;
                }
                cerr << '\n';
            }
            else{
                auto h = graph.get_handle_of_step(kmer.steps.back());
                auto length = graph.get_length(h);

                if (length > k + 1){
                    kmer.initialize(kmer.steps.back(), length - k + 1);
                    prev_has_diploid = false;
                }
            }

            prev_has_diploid = kmer.has_diploid;
        }
    });
}


int main (int argc, char* argv[]){
    path gfa_path;
    size_t k;

    CLI::App app{"App description"};

    app.add_option(
            "-i,--input_gfa",
            gfa_path,
            "Path to GFA containing phased non-overlapping segments")
            ->required();

    app.add_option(
            "-k,--kmer_size",
            k,
            "Length of kmer (k) to use")
            ->required();

    CLI11_PARSE(app, argc, argv);

    extract_haplotype_kmers_from_gfa(gfa_path, k);

    return 0;
}
