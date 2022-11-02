#include "binomial.hpp"
#include <iostream>

using std::cerr;

using gfase::binomial_coefficient;
using gfase::binomial;


int main(){

    double p = 0.5;
    for (int64_t n=0; n <= 10; n++){
        cerr << "---" << '\n';
        for (int64_t k=0; k <= n; k++){
            cerr << n << '\t' << k << '\t' << binomial_coefficient(double(n),double(k)) << '\t' << binomial(p,double(n),double(k)) << '\n';
        }
    }

    return 0;
}