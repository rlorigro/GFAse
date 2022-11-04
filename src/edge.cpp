#include "edge.hpp"

#include <algorithm>

using std::min;
using std::max;

namespace gfase {

pair<int32_t,int32_t> edge(int32_t a, int32_t b){
    if (a < b){
        return {a,b};
    }
    else{
        return {b,a};
    }
}

}