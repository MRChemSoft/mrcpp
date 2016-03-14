#ifndef TIMER_H
#define TIMER_H

#include <boost/timer/timer.hpp>

#include "TelePrompter.h"

class Timer {
public:
    Timer() : nano(1.0e9) { t.stop(); }
    void restart() { t.resume(); }
    void stop() { t.stop(); }
    double getWallTime() { return t.elapsed().wall/nano; }
    double getUserTime() { return t.elapsed().user/nano; }
    double getSystemTime() { return t.elapsed().system/nano; }
    friend std::ostream& operator<<(std::ostream &o, Timer &t) {
        int old_prec;
        GET_PRINT_PRECISION(old_prec);
        SET_PRINT_PRECISION(3);
        o << "    user   " << std::setw(6) <<  t.getUserTime();
        o << "    wall   " << std::setw(6) <<  t.getWallTime();
        SET_PRINT_PRECISION(old_prec);
        return o;
    }
private:
    const double nano;
    boost::timer::cpu_timer t;
};

#endif // TIMER_H
