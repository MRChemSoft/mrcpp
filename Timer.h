#ifndef TIMER_H
#define TIMER_H

#include <time.h>
#include <mpi.h>

#include "TelePrompter.h"

class Timer {
public:
    Timer() : time_used(0.0) { }
    Timer& operator=(const Timer &timer) {
        if (this != &timer) {
            this->time_used = timer.time_used;
        }
        return *this;
    }
    double clock_start;
    void start() {
        clock_start=gettime();
        time_used=0.0;
    }
    void restart() { clock_start=gettime();}
    void stop() { time_used += (gettime()-clock_start); }
    double getWallTime() const { return time_used; }
    double getUserTime() const { return time_used; }
    double getSystemTime() const { return time_used; }
    double gettime(){
//      return MPI_Wtime();
        timespec time1;
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
        return (double) time1.tv_sec+1.0E-9*time1.tv_nsec;
    }
    friend std::ostream& operator<<(std::ostream &o, const Timer &time_used) {
        int old_prec;
        GET_PRINT_PRECISION(old_prec);
        SET_PRINT_PRECISION(3);
        o << "    user   " << std::setw(6) <<  time_used.getUserTime();
        o << "    wall   " << std::setw(6) <<  time_used.getWallTime();
        SET_PRINT_PRECISION(old_prec);
        return o;
    }
private:
    double time_used;
};

#endif // TIMER_H
