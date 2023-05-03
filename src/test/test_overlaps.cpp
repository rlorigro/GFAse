#include "Overlaps.hpp"
#include "Filesystem.hpp"
#include "gfa_to_handle.hpp"
#include "handle_to_gfa.hpp"
#include "graph_utility.hpp"
#include "IncrementalIdMap.hpp"

#include "bdsg/hash_graph.hpp"
#include "handlegraph/util.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <cassert>
#include <sstream>
#include <algorithm>
#include <map>

using namespace std;

using ghc::filesystem::path;
using bdsg::HashGraph;
using handlegraph::handle_t;
using handlegraph::nid_t;

int ref_len(vector<pair<int, char>>& correct) {
    int l = 0;
    for (auto p : correct) {
        if (p.second != 'I') {
            l += p.first;
        }
    }
    return l;
}
int query_len(vector<pair<int, char>>& correct) {
    int l = 0;
    for (auto p : correct) {
        if (p.second != 'D') {
            l += p.first;
        }
    }
    return l;
}

vector<pair<int, char>> reverse_ops(vector<pair<int, char>>& correct) {
    auto reversed = correct;
    reverse(reversed.begin(), reversed.end());
    for (auto& p : reversed) {
        if (p.second == 'I') {
            p.second = 'D';
        }
        else if (p.second == 'D') {
            p.second = 'I';
        }
    }
    return reversed;
}

void test_cigar_internal(Cigar& cigar, vector<pair<int, char>>& correct) {
    
    bool pass = true;
    pass &= (cigar.empty() == correct.empty());
    pass &= (cigar.size() == correct.size());
    auto lens = cigar.aligned_length();
    pass &= (lens.first == ref_len(correct));
    pass &= (lens.second == query_len(correct));
    for (size_t i = 0; i < cigar.size(); ++i) {
        pass &= (cigar.at(i).type() == correct[i].second);
        pass &= (cigar.at(i).length() == correct[i].first);
    }
    
    if (cigar.empty()) {
        pass &= (cigar.get_string() == "0M");
    }
    else {
        stringstream s;
        for (auto p : correct) {
            s << p.first << p.second;
        }
        pass &= (cigar.get_string() == s.str());
    }
    
    if (!pass) {
        cerr << "cigar test fail on correct CIGAR ";
        for (auto o : correct) {
            cerr << o.first << o.second;
        }
        cerr << endl;
        cerr << "got CIGAR ";
        for (auto o : cigar) {
            cerr << o.length() << o.type();
        }
        cerr << endl;
        assert(false);
    }
}

void test_cigar(string& cigar_str, vector<pair<int, char>>& correct) {
    
    Cigar cigar(cigar_str);
    Cigar rev_cigar = cigar.reverse();
    auto rev_correct = reverse_ops(correct);
    test_cigar_internal(cigar, correct);
    test_cigar_internal(rev_cigar, rev_correct);
}

void test_gfa_overlaps_internal(HashGraph& graph, IncrementalIdMap<string>& id_map, Overlaps& overlaps,
                               vector<tuple<string, bool, string, bool, vector<pair<int, char>>>>& correct) {
    
    map<tuple<string, bool, string, bool>, vector<pair<int, char>>> answers;
    for (auto& r : correct) {
        answers[make_tuple(get<0>(r), get<1>(r), get<2>(r), get<3>(r))] = get<4>(r);
        answers[make_tuple(get<2>(r), !get<3>(r), get<0>(r), !get<1>(r))] = reverse_ops(get<4>(r));
    }
    
    assert(overlaps.is_blunt() == answers.empty());
    
    graph.for_each_edge([&](const edge_t& edge) {
        for (auto e : {edge, edge_t(graph.flip(edge.second), graph.flip(edge.first))}) {
            
            auto key = make_tuple(id_map.get_name(graph.get_id(e.first)),
                                  graph.get_is_reverse(e.first),
                                  id_map.get_name(graph.get_id(e.second)),
                                  graph.get_is_reverse(e.second));
            
            assert(overlaps.has_overlap(graph, e.first, e.second) == (bool) answers.count(key));
            
            auto cigar = overlaps.get_overlap(graph, e.first, e.second);
            if (answers.count(key)) {
                test_cigar_internal(cigar, answers.at(key));
            }
            else {
                vector<pair<int, char>> c;
                test_cigar_internal(cigar, c);
            }
        }
    });
    
    if (!correct.empty()) {
        // test remove
        auto r = correct.back();
        correct.pop_back();
        overlaps.remove_overlap(graph,
                                graph.get_handle(id_map.get_id(get<0>(r)), get<1>(r)),
                                graph.get_handle(id_map.get_id(get<2>(r)), get<3>(r)));
        
        test_gfa_overlaps_internal(graph, id_map, overlaps, correct);
    }
}


void parse_gfa(string& data_file,
               HashGraph& graph,
               IncrementalIdMap<string>& id_map,
               Overlaps& overlaps) {

    path script_path = __FILE__;
    path project_directory = script_path.parent_path().parent_path().parent_path();
    path relative_gfa_path = data_file;
    path absolute_gfa_path = project_directory / relative_gfa_path;
    
    gfa_to_handle_graph(graph, id_map, overlaps, absolute_gfa_path);
}

void test_gfa_overlaps(string& data_file,
                       vector<tuple<string, bool, string, bool, vector<pair<int, char>>>>& correct) {
    
    HashGraph graph;
    IncrementalIdMap<string> id_map;
    Overlaps overlaps;
    
    parse_gfa(data_file, graph, id_map, overlaps);
    
    test_gfa_overlaps_internal(graph, id_map, overlaps, correct);
}

int main(){

    // cigar tests
    {
        string cigar = "2D5M9I";
        vector<pair<int, char>> correct{
            {2, 'D'},
            {5, 'M'},
            {9, 'I'}
        };
        test_cigar(cigar, correct);
    }
    {
        string cigar = "5M1D3I6X1=";
        vector<pair<int, char>> correct{
            {5, 'M'},
            {1, 'D'},
            {3, 'I'},
            {6, 'X'},
            {1, '='}
        };
        test_cigar(cigar, correct);
    }
    {
        string cigar = "0M";
        vector<pair<int, char>> correct;
        test_cigar(cigar, correct);
    }
    {
        string cigar = "*";
        vector<pair<int, char>> correct;
        test_cigar(cigar, correct);
    }
    cerr << "Passed CIGAR tests!" << endl;
    
    // overlap tests
    {
        string file = "data/simple_chain.gfa";
        vector<tuple<string, bool, string, bool, vector<pair<int, char>>>> correct{
            
        };
        test_gfa_overlaps(file, correct);
    }
    {
        string file = "data/test_gfa1.gfa";
        vector<tuple<string, bool, string, bool, vector<pair<int, char>>>> correct{
            {"11", false, "12", true, {{4, 'M'}}},
            {"12", true, "13", false, {{5, 'M'}}},
            {"11", false, "13", false, {{3, 'M'}}}
        };
        test_gfa_overlaps(file, correct);
    }
    
    {
        HashGraph graph;
        IncrementalIdMap<string> id_map;
        Overlaps overlaps;
        
        string data_file = "data/test_gfa1.gfa";
        parse_gfa(data_file, graph, id_map, overlaps);\
        
        unzip(graph, id_map, overlaps);
        
        bool pass = true;
        pass &= (graph.get_node_count() == 1);
        graph.for_each_handle([&](const handle_t& h) {
            // one of the two orientations
            pass &= (graph.get_sequence(h) == "AATCAAGGT" || graph.get_sequence(h) == "ACCTTGATT");
        });
        if (!pass) {
            cerr << "unzip with stitching test fail, got:\n";
            graph.for_each_handle([&](const handle_t& h) {
                cerr << graph.get_id(h) << " " << graph.get_sequence(h) << endl;
                graph.follow_edges(h, true, [&](const handle_t& p) {
                    cerr << '\t' << graph.get_id(p) << " " << graph.get_is_reverse(p) << " <-" << endl;
                });
                graph.follow_edges(h, false, [&](const handle_t& p) {
                    cerr << "\t-> " << graph.get_id(p) << " " << graph.get_is_reverse(p) << endl;
                });
            });
            assert(false);
        }
    }
    
    cerr << "All tests passed!" << endl;
    return 0;
}
