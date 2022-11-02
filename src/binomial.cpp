#include "binomial.hpp"
#include <stdexcept>
#include <iostream>
#include <cmath>

using std::runtime_error;
using std::pow;
using std::beta;
using std::round;
using std::cerr;


namespace gfase{

double binomial_coefficient(double n, double k){
    if (k > n){
        throw runtime_error("ERROR: cannot compute binomial coefficient for k > n");
    }

    double intermediate = (1.0/(n+1.0)) * (1.0/beta(k+1.0, (n-k)+1.0));
    auto result = double(round(intermediate));

    return result;
}


double binomial(double p, double n, double k){
    auto c = binomial_coefficient(n,k);

    auto a = pow(p,k);
    auto b = pow((1-p),(n-k));

    return c*a*b;
}


}
