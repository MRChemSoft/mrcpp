#ifndef TIMER_H
#define TIMER_H

#include <boost/timer/timer.hpp>

#include "TelePrompter.h"

class Timer {
public:
    Timer() : nano(1.0e9) { t.stop(); }
    Timer& operator=(const Timer &timer) {
        if (this != &timer) {
            this->t = timer.t;
        }
        return *this;
    }
    void restart() { t.resume(); }
    void stop() { t.stop(); }
    double getWallTime() const { return t.elapsed().wall/nano; }
    double getUserTime() const { return t.elapsed().user/nano; }
    double getSystemTime() const { return t.elapsed().system/nano; }
    friend std::ostream& operator<<(std::ostream &o, const Timer &t) {
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
