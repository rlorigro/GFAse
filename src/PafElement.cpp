#include "PafElement.hpp"

namespace gfase{

PafElement::PafElement(string& target_name,
                       string& query_name,
                       uint32_t start,
                       uint32_t stop,
                       uint32_t map_quality,
                       bool is_reverse):
    target_name(target_name),
    query_name(query_name),
    start(start),
    stop(stop),
    map_quality(map_quality),
    is_reverse(is_reverse)
{}

ostream& operator<<(ostream& o, PafElement& e){
    cerr << e.query_name << " "
         << e.is_reverse << " "
         << e.target_name << " "
         << e.start << " "
         << e.stop << " "
         << e.map_quality << " ";

    return o;
}

}