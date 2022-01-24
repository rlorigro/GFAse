#include "Color.hpp"
#include <iostream>
using std::cerr;


int main(){
    double r = 0;
    double g = 0.5;
    double b = 1;

    cerr << rgb_to_hex(r,g,b) << '\n';

    return 0;
}
