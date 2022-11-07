#include "Chainer.hpp"

#include <queue>

using std::queue;


namespace gfase{


void Chainer::get_chain(
        const HandleGraph& graph,
        const nid_t& start_node,
        deque <set <nid_t> >& chain
        ){

    unordered_set<nid_t> visited;
    queue <set <nid_t> > q;

    // Initialize start node as an item in the chain
    // Only chainable nodes should be considered
    {
        bool is_diploid = diploid_nodes.count(start_node);
        bool is_haploid = haploid_nodes.count(start_node);

        // Check that this has NOT been visited before queuing it
        if (is_diploid or is_haploid) {
            q.emplace();
            q.back().emplace(start_node);

            // If this node is diploid, also add the alt/pair but don't bother to queue it.
            if (is_diploid){
                auto alt_node = node_pairs.at(start_node);
                visited.emplace(alt_node);
                q.back().emplace(alt_node);
            }
        }
    }

    // Search left
    while (not q.empty()) {
        set <nid_t> item = q.front();
        q.pop();

        chain.emplace_front(item);

        auto h = graph.get_handle(*item.begin());

        graph.follow_edges(h, true, [&](const handle_t& other_handle) {
            auto next_node = graph.get_id(other_handle);

            // Attempt to add this to the set of visited nodes
            auto unvisited = visited.emplace(next_node).second;

            // Only chainable nodes should be considered
            bool is_diploid = diploid_nodes.count(next_node);
            bool is_haploid = haploid_nodes.count(next_node);

            // Check that this has NOT been visited before queuing it
            if (unvisited and (is_diploid or is_haploid)) {
                q.emplace();
                q.back().emplace(next_node);

                // If this node is diploid, also add the alt/pair but don't bother to queue it.
                if (is_diploid){
                    auto alt_node = node_pairs.at(next_node);
                    visited.emplace(alt_node);
                    q.back().emplace(alt_node);
                }
            }
        });
    }

    // Initialize start node as an item in the chain
    // Only chainable nodes should be considered
    {
        bool is_diploid = diploid_nodes.count(start_node);
        bool is_haploid = haploid_nodes.count(start_node);

        // Check that this has NOT been visited before queuing it
        if (is_diploid or is_haploid) {
            q.emplace();
            q.back().emplace(start_node);

            // If this node is diploid, also add the alt/pair but don't bother to queue it.
            if (is_diploid){
                auto alt_node = node_pairs.at(start_node);
                visited.emplace(alt_node);
                q.back().emplace(alt_node);
            }
        }
    }

    // Search right
    while (not q.empty()) {
        set <nid_t> item = q.front();
        q.pop();

        if (item.count(start_node) == 0) {
            chain.emplace_back(item);
        }

        auto h = graph.get_handle(*item.begin());

        graph.follow_edges(h, false, [&](const handle_t& other_handle) {
            auto next_node = graph.get_id(other_handle);

            // Attempt to add this to the set of visited nodes
            auto unvisited = visited.emplace(next_node).second;

            // Only chainable nodes should be considered
            bool is_diploid = diploid_nodes.count(next_node);
            bool is_haploid = haploid_nodes.count(next_node);

            // Check that this has NOT been visited before queuing it
            if (unvisited and (is_diploid or is_haploid)) {
                q.emplace();
                q.back().emplace(next_node);

                // If this node is diploid, also add the alt/pair but don't bother to queue it.
                if (is_diploid){
                    auto alt_node = node_pairs.at(next_node);
                    visited.emplace(alt_node);
                    q.back().emplace(alt_node);
                }
            }
        });
    }
}


void Chainer::for_each_chain(
        HandleGraph& graph,
        const function<void(chain_t& chain)>& f
        ){

    unordered_set<nid_t> visited;

    graph.for_each_handle([&](const handle_t& h){
        auto id = graph.get_id(h);

        if (visited.count(id)){
            return;
        }

        chain_t chain;

        get_chain(graph, id, chain);

        // Mark all the nodes in the chain as visited
        for (auto& item: chain){
            for (auto& n: item){
                visited.emplace(n);
            }
        }

        if (not chain.empty()) {
            f(chain);
        }
    });
}


size_t Chainer::get_path_length(const PathHandleGraph& graph, const path_handle_t& p) const{
    size_t length = 0;
    graph.for_each_step_in_path(p, [&](const step_handle_t& s) {
        length++;
    });

    return length;
}


void Chainer::new_paths(size_t chain_index, MutablePathDeletableHandleGraph& graph, array<path_handle_t,2>& paths, bool check_empty){
    if (check_empty) {
        size_t total_steps = 0;
        for (auto& p: paths) {
            total_steps += get_path_length(graph, p);
        }

        // Don't make new pair of paths if previous one was unused
        if (total_steps == 0) {
            return;
        }
    }

    string name_0 = "gfase_hap_" + to_string(chain_index) +"_0";
    string name_1 = "gfase_hap_" + to_string(chain_index) +"_1";

    cerr << name_0 << ' ' << name_1 << '\n';

    paths[0] = graph.create_path_handle(name_0);
    paths[1] = graph.create_path_handle(name_1);
}


bool Chainer::process_diploid_chain_element(
        const set<nid_t>& chain_element,
        array<path_handle_t,2>& paths,
        MutablePathDeletableHandleGraph& graph,
        const MultiContactGraph& contact_graph
        ){

    auto a = *chain_element.begin();
    auto b = *chain_element.rbegin();
    array<nid_t,2> nodes = {-1,-1};

    bool phasable = false;

    if (contact_graph.has_node(int32_t(a)) and contact_graph.has_node(int32_t(b))){
        auto phase_a = contact_graph.get_partition(int32_t(a));
        auto phase_b = contact_graph.get_partition(int32_t(b));

        if (phase_a != phase_b and phase_a != 0 and phase_b != 0) {
            // from  -1 or 1
            // to     0 or 1
            nodes[phase_a == 1] = a;
            nodes[phase_b == 1] = b;

            graph.append_step(paths[0], graph.get_handle(nodes[0]));
            graph.append_step(paths[1], graph.get_handle(nodes[1]));

            phasable = true;
        }
    }

    return phasable;
}


void Chainer::process_haploid_chain_element(
        const set<nid_t>& chain_element,
        array<path_handle_t,2>& paths,
        MutablePathDeletableHandleGraph& graph
        ){
    auto n = *chain_element.begin();

    graph.append_step(paths[0], graph.get_handle(n));
    graph.append_step(paths[1], graph.get_handle(n));
}


void Chainer::generate_chain_paths(
        MutablePathDeletableHandleGraph& graph,
        const IncrementalIdMap<string>& id_map,
        const MultiContactGraph& contact_graph
        ){

    find_chainable_nodes(graph, id_map);

    size_t c = 0;
    for_each_chain(graph, [&](chain_t& chain){
        // Singletons don't need to be unzipped
        if (chain.size() == 1){
            return;
        }

        array<path_handle_t,2> paths;
        new_paths(c, graph, paths, false);

        for (auto& item: chain){
            cerr << c << ' ';
            for (auto n: item){
                cerr << id_map.get_name(n) << ' ';
            }
            cerr << '\n';

            if (item.size() == 2){
                bool phasable = process_diploid_chain_element(item, paths, graph, contact_graph);

                // If unphasable, don't append unannotated nodes to the path, don't try to assign phase,
                // start new path.
                if (not phasable) {
                    cerr << "Warning: breaking chain at unphasable item" << '\n';
                    new_paths(c, graph, paths, true);
                    c++;
                }
            }
            else if (item.size() == 1){
                process_haploid_chain_element(item, paths, graph);
            }
            else{
                for (auto& id: item){
                    cerr << id_map.get_name(id) << '\n';
                }
                throw runtime_error("ERROR: unchainable items in chain, see above for details");
            }
        }

        c++;
    });

    vector <path_handle_t> to_be_destroyed;

    graph.for_each_path_handle([&](const path_handle_t& p){
        auto length = get_path_length(graph,p);

        if (length < 2){
            to_be_destroyed.emplace_back(p);
        }

        cerr << graph.get_path_name(p) << ' ';
        graph.for_each_step_in_path(p, [&](const step_handle_t& s){
            auto id = graph.get_id(graph.get_handle_of_step(s));
            auto name = id_map.get_name(id);
            cerr << name << ' ';
        });
        cerr << '\n';
    });

    for (auto& p: to_be_destroyed){
        graph.destroy_path(p);
    }
}


void Chainer::find_chainable_nodes(const HandleGraph& graph, const IncrementalIdMap<string>& id_map){
    graph.for_each_handle([&](const handle_t& h0) {
        auto id0 = graph.get_id(h0);

        // Do a two-edge walk right/left and left/right
        unordered_set<nid_t> left_second_degree_neighbors;
        unordered_set<nid_t> right_second_degree_neighbors;

        unordered_set<nid_t> left_first_degree_neighbors;
        unordered_set<nid_t> right_first_degree_neighbors;

        bool has_self_edge = false;

        graph.follow_edges(h0, true, [&](const handle_t& h1){
            auto id1 = graph.get_id(h1);
            if (id0 != id1) {
                left_first_degree_neighbors.emplace(id1);
            }
            else{
                has_self_edge = true;
            }

            graph.follow_edges(h1, false, [&](const handle_t& h2){
                auto id2 = graph.get_id(h2);

                if (id0 != id2) {
                    left_second_degree_neighbors.emplace(id2);
                }
            });
        });

        graph.follow_edges(h0, false, [&](const handle_t& h1){
            auto id1 = graph.get_id(h1);
            if (id0 != id1) {
                right_first_degree_neighbors.emplace(id1);
            }
            else{
                has_self_edge = true;
            }

            graph.follow_edges(h1, true, [&](const handle_t& h2){
                auto id2 = graph.get_id(h2);

                if (id0 != id2) {
                    right_second_degree_neighbors.emplace(id2);
                }
            });
        });

        bool is_haploid = (left_second_degree_neighbors.empty() and right_second_degree_neighbors.empty());

        bool is_symmetrical_bubble = (right_second_degree_neighbors == left_second_degree_neighbors);
        bool is_diploid_bubble = (right_second_degree_neighbors.size() == 1) and (left_second_degree_neighbors.size() == 1);

        bool is_chainable = (left_first_degree_neighbors.size() < 3) and (right_first_degree_neighbors.size() < 3);

        bool is_tip = ((left_first_degree_neighbors.empty() and (right_second_degree_neighbors.size() == 1))
                   or (right_first_degree_neighbors.empty() and (left_second_degree_neighbors.size() == 1)));


        if (is_haploid and not has_self_edge){
            haploid_nodes.emplace(graph.get_id(h0));
        }
        else if (is_symmetrical_bubble and is_diploid_bubble and is_chainable and not has_self_edge){
            nid_t id_a = graph.get_id(h0);
            nid_t id_b = *left_second_degree_neighbors.begin();

            diploid_nodes.emplace(id_a);
            diploid_nodes.emplace(id_b);
            node_pairs.emplace(id_a,id_b);
            node_pairs.emplace(id_b,id_a);
        }
        else if (is_tip and is_chainable and not has_self_edge){
            nid_t id_a = graph.get_id(h0);
            diploid_tip_nodes.emplace(id_a);
        }

//        cerr << id_map.get_name(graph.get_id(h0)) << '\n';
//        cerr << "is_symmetrical_bubble: " << is_symmetrical_bubble << '\n';
//        cerr << "is_diploid_bubble: " << is_diploid_bubble << '\n';
//        cerr << "is_chainable: " << is_chainable << '\n';
//        cerr << "is_tip: " << is_tip << '\n';
//        cerr << "is_haploid: " << is_haploid << '\n' << '\n';
    });
}


void Chainer::write_chainable_nodes_to_bandage_csv(path output_dir, const IncrementalIdMap<string>& id_map) const{
    path output_path = output_dir / "chainable_nodes.csv";
    ofstream file(output_path);

    if (not file.is_open() or not file.good()) {
        throw std::runtime_error("ERROR: could not write to file: " + output_path.string());
    }

    file << "Name" << ',' << "Color" << '\n';
    for (auto& item: haploid_nodes){
        file << id_map.get_name(item) << ',' << "Cornflower Blue" << '\n';
    }
    for (auto& item: diploid_nodes){
        file << id_map.get_name(item) << ',' << "Midnight Blue" << '\n';
    }
    for (auto& item: diploid_tip_nodes){
        file << id_map.get_name(item) << ',' << "Slate Blue" << '\n';
    }
}


}
