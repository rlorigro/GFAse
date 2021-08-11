#ifndef GFASE_INCREMENTALIDMAP_HPP
#define GFASE_INCREMENTALIDMAP_HPP

#include "handlegraph/handle_graph.hpp"
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <memory>

using std::string;
using std::vector;
using std::unordered_map;
using handlegraph::nid_t;
using std::unique_ptr;
using std::runtime_error;
using std::unique_ptr;
using std::shared_ptr;
using std::make_unique;
using std::make_shared;
using std::to_string;


namespace gfase{

template <class T> class IncrementalIdMap {
public:
    /// Attributes ///

    vector <T> names;
    unordered_map <T, int64_t> ids;
    bool zero_based;

    /// Methods ///

    IncrementalIdMap<T>(bool zero_based=false);

    // Add a node ID to the running list, do whatever needs to be done to make sure the mapping is reversible, and then
    // return its incremental ID, based on the number of nodes added so far
    int64_t insert(const T& s);

    // Find the original node ID from its integer ID
    T get_name(int64_t id) const;
    int64_t get_id(const T& name) const;

    // Check if key/value has been added already, returns true if it exists
    bool exists(const T& name);
    bool exists(int64_t id);
};


template<class T> IncrementalIdMap<T>::IncrementalIdMap(bool zero_based):
    zero_based(zero_based)
{}


template <class T> int64_t IncrementalIdMap<T>::insert(const T& s) {
    if (exists(s)){
        throw runtime_error("Error: attempted to insert duplicate key with id: " + to_string(get_id(s)));
    }

    // Make a copy of the node name string
    names.push_back(s);

    // Create an integer node ID (starting from 1)
    int64_t id = names.size() - zero_based;

    // Create a reverse mapping
    ids.emplace(s,id);

    // For convenience, return the ID number that was generated
    return id;
}


template<class T> bool IncrementalIdMap<T>::exists(const T& name){
    return (ids.find(name) != ids.end());
}


template<class T> bool IncrementalIdMap<T>::exists(int64_t id){
    return (id >= 0 and id <= names.size() - zero_based);
}


template<class T> int64_t IncrementalIdMap<T>::get_id(const T& name) const{
    return ids.at(name);
}


template<class T> T IncrementalIdMap<T>::get_name(int64_t id) const{
    return names[id-1+zero_based];
}


}
#endif //GFASE_INCREMENTALID_HPP
