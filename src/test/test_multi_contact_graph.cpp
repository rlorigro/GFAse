#include "MultiContactGraph.hpp"

using gfase::NonBipartiteEdgeException;
using gfase::MultiContactGraph;
using gfase::alt_component_t;

#include <iostream>

using std::exception;
using std::cerr;


int main(){
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

        cerr << "same side 0 vs 1: " << int(g.of_same_component_side(0,1)) << '\n';
        cerr << "same side 0 vs 2: " << int(g.of_same_component_side(0,2)) << '\n';
        cerr << "same side 0 vs 3: " << int(g.of_same_component_side(0,3)) << '\n';
        cerr << "same side 0 vs 4: " << int(g.of_same_component_side(0,4)) << '\n';
    }
    cerr << "TESTING triangle component:" << '\n';
    {
        MultiContactGraph g;

        g.insert_node(0);
        g.insert_node(1);
        g.insert_node(2);

        g.add_alt(0,1);
        g.add_alt(1,2);

        try {
            g.add_alt(2,0);
        }
        catch(NonBipartiteEdgeException& e){
            cerr << e.what() << '\n';
            cerr << "successfully caught invalid component" <<'\n';

            cerr << "Testing exception data: " << '\n';

            for (auto& item: e.component_a.first){
                cerr << "0 " << item << '\n';
            }
            for (auto& item: e.component_a.second){
                cerr << "1 " << item << '\n';
            }

            for (auto& item: e.component_b.first){
                cerr << "0 " << item << '\n';
            }
            for (auto& item: e.component_b.second){
                cerr << "1 " << item << '\n';
            }
            cerr << '\n';
        }

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
        catch (NonBipartiteEdgeException& e){
            cerr << e.what() << '\n';
            cerr << "successfully caught invalid component" << '\n';
            cerr << "Testing exception data: " << '\n';

            for (auto& item: e.component_a.first){
                cerr << "0 " << item << '\n';
            }
            for (auto& item: e.component_a.second){
                cerr << "1 " << item << '\n';
            }

            for (auto& item: e.component_b.first){
                cerr << "0 " << item << '\n';
            }
            for (auto& item: e.component_b.second){
                cerr << "1 " << item << '\n';
            }
        }

        g.validate_alts();
    }

    cerr << "TESTING component merge:" << '\n';
    {
        MultiContactGraph g;

        g.insert_node(0);
        g.set_partition(0,-1);
        g.insert_node(1);
        g.set_partition(1,1);
        g.insert_node(2);
        g.set_partition(2,-1);
        g.insert_node(3);
        g.set_partition(3,1);
        g.insert_node(4);
        g.set_partition(4,-1);
        g.insert_node(5);
        g.set_partition(5,1);

        g.add_alt(0,1);
        g.add_alt(2,3);
        g.add_alt(4,5);

        alt_component_t a;
        alt_component_t b;
        alt_component_t c;

        g.get_alt_component(0, false, a);
        g.get_alt_component(2, false, b);

        cerr << "a" << '\n';
        for (auto& item: a.first){
            cerr << '\t' << 0 << " | " << item << ' ' << int(g.get_partition(item)) << '\n';
        }
        for (auto& item: a.second){
            cerr << '\t' << 1 << " | " << item << ' ' << int(g.get_partition(item)) << '\n';
        }
        cerr << "b" << '\n';
        for (auto& item: b.first){
            cerr << '\t' << 0 << " | " << item << ' ' << int(g.get_partition(item)) << '\n';
        }
        for (auto& item: b.second){
            cerr << '\t' << 1 << " | " << item << ' ' << int(g.get_partition(item)) << '\n';
        }

        g.add_alt(a,b);
        g.get_alt_component(0,false,c);

        cerr << "c" << '\n';
        for (auto& item: c.first){
            cerr << '\t' << 0 << " | " << item << ' ' << int(g.get_partition(item)) << '\n';
        }
        for (auto& item: c.second){
            cerr << '\t' << 1 << " | " << item << ' ' << int(g.get_partition(item)) << '\n';
        }

        cerr << "same side 0 vs 1: " << int(g.of_same_component_side(0,1)) << '\n';
        cerr << "same side 0 vs 2: " << int(g.of_same_component_side(0,2)) << '\n';
        cerr << "same side 0 vs 3: " << int(g.of_same_component_side(0,3)) << '\n';
        cerr << "same side 0 vs 4: " << int(g.of_same_component_side(0,4)) << '\n';


    }


    return 0;
}
