#include "Timer.hpp"
#include <iostream>
#include <thread>         // std::this_thread::sleep_for

using gfase::Timer;
using std::this_thread::sleep_for;
using std::cerr;


int main() {
    Timer t;

    cerr << t << "start" << '\n';

    auto d = std::chrono::seconds(5);
    std::this_thread::sleep_for(d);

    cerr << t << "finish" << '\n';

    return 0;
}

