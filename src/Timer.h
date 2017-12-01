#pragma once

#include <chrono>
#include <iostream>

typedef std::chrono::time_point<std::chrono::high_resolution_clock> timeT;

class Timer {
public:
    Timer(bool start_timer = true);
    Timer& operator=(const Timer &timer);

    void start();
    void resume();
    void stop();

    double getWallTime() const;

    friend std::ostream& operator<<(std::ostream &o, const Timer &timer);
private:
    bool running;
    double time_used;
    timeT clock_start;

    timeT now();
    double diffTime(timeT t2, timeT t1);
};
