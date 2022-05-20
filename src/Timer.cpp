#include "Timer.hpp"

using std::to_string;

namespace gfase{

Timer::Timer():
        start(std::chrono::steady_clock::now())
{}


string Timer::elapsed() const {
    auto d = std::chrono::steady_clock::now() - start;

    const auto h = duration_cast<hours>(d);
    const auto m = duration_cast<minutes>(d - h);
    const auto s = duration_cast<seconds>(d - h - m);

    string ds;
    ds.append("[");
    ds.append(to_string(h.count()));
    ds.append("h ");
    ds.append(to_string(m.count()));
    ds.append("m ");
    ds.append(to_string(s.count()));
    ds.append("s");
    ds.append("] ");

    return ds;
}


void Timer::reset() {
    start = std::chrono::steady_clock::now();
}


}

ostream& operator<<(ostream& o, const gfase::Timer& t) {
    o << t.elapsed();

    return o;
}


