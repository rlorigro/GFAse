#include "HaplotypePathKmer.hpp"

#include <stdexcept>
#include <string>

using std::runtime_error;
using std::to_string;
using std::string;
using std::cerr;


namespace gfase{

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
        auto p = graph.get_path_handle_of_step(s);

        // If this is the terminal step, its ok for it to overlap (because the paths may have been extended)
        // If internal node is multiploid, then it must be error
        if (not (s == graph.path_back(p) or s == graph.path_begin(p))) {
//            auto path_name = graph.get_path_name(p);
//
//            cerr << "ERROR: node has >2 paths: " << '\n';
//
//            graph.for_each_step_on_handle(h, [&](const step_handle_t& s_other){
//                cerr << '\t' << graph.get_path_name(graph.get_path_handle_of_step(s_other)) << '\n';
//            });
//
//            throw runtime_error(
//                    "ERROR: multiploid step in path " + path_name + " for node: " + to_string(graph.get_id(h)));
        }
        else{
            // If terminal node and multiploid, it should not be treated as diploid
            result = false;
        }
    }

    return result;
}


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

    step_is_diploid = is_haplotype_bubble(graph, s);

    steps.emplace_back(s);
    lengths.emplace_back(graph.get_length(h));
    is_diploid.emplace_back(step_is_diploid);

    if (lengths.back() > 0) {
        sequence.emplace_back(graph.get_base(h, index));
    }

    update_has_diploid();

    start_index = index;
    stop_index = index;

    terminal_step = graph.path_back(path);
    terminal_index = graph.get_length(graph.get_handle_of_step(terminal_step));

    while (sequence.size() < k - 1){
        // Preload the kmer queue up until it almost has k elements
        this->step();

        // Verify that this kmer actually fits in the path
        if (steps.back() == terminal_step and stop_index + 1 >= terminal_index) {
            // Insufficient path length, kmer will be shorter than k
            break;
        }
    }

//    if (not sufficient_path_length){
//        throw runtime_error("ERROR: path " + graph.get_path_name(path) + " does not have sufficient length for kmer size: " + to_string(k));
//    }

}


void HaplotypePathKmer::for_each_haploid_kmer(const function<void(deque<char>& sequence)>& f){
    while (step()) {
        if (has_diploid) {
            f(sequence);
        }
        else {
            auto h = graph.get_handle_of_step(steps.back());
            auto length = graph.get_length(h);

            // Skip to the end of long nodes if they are not diploid
            if (length > k + 1) {
                initialize(steps.back(), length - k + 1);
            }
        }
    }
}


void HaplotypePathKmer::print(){
    cerr << "start: " << start_index << '\n';
    cerr << "stop: " << stop_index << '\n';

    cerr << "L: ";
    for (auto& l: lengths){
        cerr << l << ',';
    }
    cerr << '\n';

    cerr << "S: ";
    for (auto& s: steps){
        cerr << graph.get_id(graph.get_handle_of_step(s)) << ',';
    }
    cerr << '\n';

    cerr << "D: ";
    for (auto& d: is_diploid){
        cerr << d << ',';
    }
    cerr << '\n';

    cerr << "K: ";
    for (auto& c: sequence){
        cerr << string(1,c) << ',';
    }
    cerr << '\n';
    cerr << '\n';

}


// TODO: switch to standard queue (only need to pop front and push back)
bool HaplotypePathKmer::step(){
    bool has_next_step = true;
    bool found_empty_node = false;
    bool moved = false;

    // If it is safe to increment this node
    if (stop_index + 1 < lengths.back()){
        stop_index++;

        auto h = graph.get_handle_of_step(steps.back());
        sequence.emplace_back(graph.get_base(h,stop_index));

        // Use the deque like a cyclic queue
        if (sequence.size() > k) {
            sequence.pop_front();
            moved = true;
        }
    }
    // If the iterator is at the end of the node
    else {
        if (steps.back() != terminal_step){
            step_handle_t next_step = graph.get_next_step(steps.back());
            auto next_handle = graph.get_handle_of_step(next_step);

            bool step_is_diploid = is_haplotype_bubble(graph, next_step);

            steps.emplace_back(next_step);
            lengths.emplace_back(graph.get_length(next_handle));
            is_diploid.emplace_back(step_is_diploid);

            if (lengths.back() > 0) {
                sequence.emplace_back(graph.get_base(next_handle, 0));
            }
            else{
                found_empty_node = true;
            }

            update_has_diploid();

            // Use the deque like a cyclic queue
            if (sequence.size() > k) {
                sequence.pop_front();
                moved = true;
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
    if (start_index + 1 < lengths.front()) {
        if (moved) {
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

    if (found_empty_node){
        return step();
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


step_handle_t HaplotypePathKmer::get_step_of_kmer_start() const{
    return steps.front();
}


step_handle_t HaplotypePathKmer::get_step_of_kmer_end() const{
    return steps.back();
}


size_t HaplotypePathKmer::get_index_of_kmer_start() const{
    return start_index;
}


}