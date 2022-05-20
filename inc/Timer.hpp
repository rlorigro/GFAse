#ifndef GFASE_DURATION_HPP
#define GFASE_DURATION_HPP

#include <iostream>
#include <chrono>

using std::chrono::duration_cast;
using std::chrono::steady_clock;
using std::chrono::duration;

using std::chrono::hours;
using std::chrono::minutes;
using std::chrono::seconds;
using std::chrono::milliseconds;

using std::ostream;
using std::string;


namespace gfase{

class Timer{
    steady_clock::time_point start;

public:
    Timer();
    string elapsed() const;
    void reset();
};

}

ostream& operator<<(ostream& o, const gfase::Timer& t);


#endif //GFASE_DURATION_HPP
