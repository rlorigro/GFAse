#ifndef GFASE_PAFELEMENT_HPP
#define GFASE_PAFELEMENT_HPP

#include <string>
#include <iostream>

using std::string;
using std::ostream;
using std::cerr;

namespace gfase {


class PafElement {
public:
    string target_name;
    string query_name;
    uint32_t start;
    uint32_t stop;
    uint32_t map_quality;
    bool is_reverse;

    PafElement()=default;
    PafElement(string& target_name,
               string& query_name,
               uint32_t start,
               uint32_t stop,
               uint32_t map_quality,
               bool is_reverse);
};

ostream& operator<<(ostream& o, PafElement& e);

}

#endif //GFASE_PAFELEMENT_HPP
