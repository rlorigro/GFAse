#include "OverlapMap.hpp"

namespace gfase {

void RegionalOverlapMap::insert(string& region_name, uint32_t start, uint32_t stop, uint32_t id) {
    // Check if this region has been encountered before, and initialize it if needed
    if (intervals.count(region_name) == 0) {
        interval_map<uint32_t, set<uint32_t>, total_enricher> empty_interval_map;
        intervals.emplace(region_name, empty_interval_map);
    }

    // Construct a Boost interval for this numeric range/coord
    auto a = interval<uint32_t>::right_open(start, stop);
    set<uint32_t> b = {id};

    // Within this contig, add the numeric interval
    intervals.at(region_name).add({a, b});
    size++;
}


void RegionalOverlapMap::print(ostream& out) {
    for (auto& item: intervals) {
        out << item.first << '\n';
        for (auto& item1: item.second) {
            out << item1.first << " -> ";
            for (auto& item2: item1.second) {
                out << item2 << " ";
            }
            out << '\n';
        }
    }
}


RegionalOverlapMap::RegionalOverlapMap() :
        intervals(),
        size(0)
{}

}