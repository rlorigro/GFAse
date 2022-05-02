#ifndef GFASE_OVERLAPMAP_HPP
#define GFASE_OVERLAPMAP_HPP

#include "boost/icl/split_interval_map.hpp"
#include "boost/icl/interval_map.hpp"
#include "boost/icl/interval.hpp"

using boost::icl::partial_enricher;
using boost::icl::total_enricher;
using boost::icl::partial_absorber;
using boost::icl::total_absorber;
using boost::icl::split_interval_map;
using boost::icl::interval_map;
using boost::icl::interval;

#include <unordered_map>
#include <ostream>
#include <string>
#include <set>

using std::unordered_map;
using std::ostream;
using std::string;
using std::set;

namespace gfase {

class RegionalOverlapMap {
public:
    unordered_map<string, interval_map<uint32_t, set<uint32_t>, total_enricher> > intervals;
    size_t size;

    void insert(string& region_name, uint32_t start, uint32_t stop, uint32_t id);

    RegionalOverlapMap();

    void print(ostream& out);
};

}

#endif //GFASE_OVERLAPMAP_HPP
