#ifndef GFASE_INCREMENTALIDMAP_HPP
#define GFASE_INCREMENTALIDMAP_HPP

#include "Filesystem.hpp"

using ghc::filesystem::path;

#include <ostream>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <memory>

using std::ofstream;
using std::ifstream;
using std::string;
using std::vector;
using std::unordered_map;
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
    IncrementalIdMap<T>(path csv_path);

    // Add a node ID to the running list, do whatever needs to be done to make sure the mapping is reversible, and then
    // return its incremental ID, based on the number of nodes added so far
    int64_t insert(const T& s);
    int64_t try_insert(const T& s);
    int64_t do_insert(const T& s);

    // Find the original node ID from its integer ID
    T get_name(int64_t id) const;
    int64_t get_id(const T& name) const;

    // Check if key/value has been added already, returns true if it exists
    bool exists(const T& name) const;
    bool exists(int64_t id) const;

    void write_to_csv(path output_path) const;
    size_t size() const;
};


template<class T> void IncrementalIdMap<T>::write_to_csv(path output_path) const{
    ofstream file(output_path);

    if (not (file.is_open() and file.good())){
        throw runtime_error("ERROR: could not write to file: " + output_path.string());
    }

    for (const auto& name: names){
        file << get_id(name) << ',' << name << '\n';
    }
}


template<class T> IncrementalIdMap<T>::IncrementalIdMap(path csv_path){
    ifstream file(csv_path);

    if (not (file.is_open() and file.good())){
        throw runtime_error("ERROR: could not read file: " + csv_path.string());
    }

    char c;
    string name;
    string id_token;
    size_t n_delimiters = 0;
    size_t n_lines = 0;

    while (file.get(c)){
        if (c == '\n'){
            if (n_lines == 0){
                size_t id = stoll(id_token);
                if (id == 0){
                    zero_based = true;
                }
                else if (id == 1){
                    zero_based = false;
                }
                else{
                    throw runtime_error("ERROR: first id is not 0 or 1: " + csv_path.string());
                }
            }

            insert(name);
            n_lines++;

            n_delimiters = 0;

            name.clear();
            id_token.clear();
        }
        else if (c == ','){
            n_delimiters++;
        }
        else{
            if (n_delimiters == 0){
                id_token += c;
            }
            else if (n_delimiters == 1){
                name += c;
            }
            else{
                throw runtime_error("ERROR: too many delimiters for line in file: " + csv_path.string());
            }
        }
    }
}


template<class T> size_t IncrementalIdMap<T>::size() const{
    return names.size();
}


template<class T> IncrementalIdMap<T>::IncrementalIdMap(bool zero_based):
    zero_based(zero_based)
{}


template <class T> int64_t IncrementalIdMap<T>::insert(const T& s) {
    if (exists(s)){
        throw runtime_error("Error: attempted to insert duplicate key with id: " + to_string(get_id(s)));
    }

    return do_insert(s);
}


template <class T> int64_t IncrementalIdMap<T>::try_insert(const T& s) {
    auto result = ids.find(s);

    if (result == ids.end()){
        return do_insert(s);
    }
    else{
        return result->second;
    }
}


template <class T> int64_t IncrementalIdMap<T>::do_insert(const T& s) {
    // Make a copy of the node name string
    names.push_back(s);

    // Create an integer node ID (starting from 1)
    int64_t id = names.size() - zero_based;

    // Create a reverse mapping
    ids.emplace(s,id);

    // For convenience, return the ID number that was generated
    return id;
}




template<class T> bool IncrementalIdMap<T>::exists(const T& name) const{
    return (ids.find(name) != ids.end());
}


template<class T> bool IncrementalIdMap<T>::exists(int64_t id) const{
    return (id >= 0 and id <= names.size() - zero_based);
}


template<class T> int64_t IncrementalIdMap<T>::get_id(const T& name) const{
    return ids.at(name);
}


template<class T> T IncrementalIdMap<T>::get_name(int64_t id) const{
    return names.at(id-1+zero_based);
}


}
#endif //GFASE_INCREMENTALID_HPP
