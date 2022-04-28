#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <set>


using std::set_intersection;
using std::string;
using std::vector;
using std::cerr;
using std::set;


int main(){
    {
        set<string> a = {"a", "b", "c"};
        set<string> b = {"b", "c", "d"};

        vector<string> result;

        set_intersection(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(result));

        for (auto& item: result) {
            cerr << item << '\n';
        }
    }

    {
        set<string> a = {};
        set<string> b = {"b", "c", "d"};

        vector<string> result;

        set_intersection(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(result));

        for (auto& item: result) {
            cerr << item << '\n';
        }
    }

    return 0;
}