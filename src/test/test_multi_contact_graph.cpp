#include "MultiContactGraph.hpp"

using gfase::MultiContactGraph;
using gfase::alt_component_t;

#include <iostream>

using std::exception;
using std::cerr;


int main(){
    cerr << "TESTING triangle component:" << '\n';
    {
        MultiContactGraph g;

        g.insert_node(0);
        g.insert_node(1);
        g.insert_node(2);

        g.add_alt(0,1);
        g.add_alt(1,2);
        g.add_alt(2,0);

        alt_component_t c;

        g.get_alt_component(0, false, c);

        for (auto& item: c.first){
            cerr << 0 << ' ' << item << '\n';
        }
        for (auto& item: c.second){
            cerr << 1 << ' ' << item << '\n';
        }

        try {
            g.get_alt_component(0, true, c);
        }
        catch (exception& e){
            cerr << e.what() << '\n';
            cerr << "successfully caught invalid component" <<'\n';
        }
    }
    cerr << "TESTING normal component:" << '\n';
    {
        MultiContactGraph g;

        g.insert_node(0);
        g.insert_node(1);
        g.insert_node(2);

        g.add_alt(0,1);
        g.add_alt(0,2);

        alt_component_t c;

        g.get_alt_component(0, false, c);

        for (auto& item: c.first){
            cerr << 0 << ' ' << item << '\n';
        }
        for (auto& item: c.second){
            cerr << 1 << ' ' << item << '\n';
        }
    }
    cerr << "TESTING zigzag component:" << '\n';
    {
        MultiContactGraph g;

        g.insert_node(0);
        g.insert_node(1);
        g.insert_node(2);
        g.insert_node(3);
        g.insert_node(4);

        g.add_alt(0,1);
        g.add_alt(1,2);
        g.add_alt(2,3);
        g.add_alt(3,4);

        alt_component_t c;

        g.get_alt_component(0, false, c);

        for (auto& item: c.first){
            cerr << 0 << ' ' << item << '\n';
        }
        for (auto& item: c.second){
            cerr << 1 << ' ' << item << '\n';
        }
    }

    return 0;
}
