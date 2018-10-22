#pragma once

#include <chrono>
#include <iostream>

namespace mrcpp {

typedef std::chrono::time_point<std::chrono::high_resolution_clock> timeT;

class Timer final {
public:
    Timer(bool start_timer = true);
    Timer(const Timer &timer);
    Timer& operator=(const Timer &timer);

    void start();
    void resume();
    void stop();

    double getWallTime() const;

    friend std::ostream& operator<<(std::ostream &o, const Timer &timer) { return timer.print(o); }
private:
    bool running;
    double time_used;
    timeT clock_start;

    timeT now();
    double diffTime(timeT t2, timeT t1);

    std::ostream& print(std::ostream &o) const;
};

}
