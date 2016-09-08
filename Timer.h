#ifndef TIMER_H
#define TIMER_H

//#include <boost/timer/timer.hpp>
#include <time.h>
#include <chrono>
#include <mpi.h>

#include "TelePrompter.h"

class Timer {
public:
  //    Timer() : nano(1.0e9) { t.stop(); }
 Timer() {start(); }
    Timer& operator=(const Timer &timer) {
        if (this != &timer) {
            this->time_used = timer.time_used;
        }
        return *this;
    }
    //    double clock_start;
    //std::chrono::system_clock::time_point clock_start = std::chrono::system_clock::now();
    /*    void restart() { t.resume(); }
    void stop() { t.stop(); }
    double getWallTime() const { return t.elapsed().wall/nano; }
    double getUserTime() const { return t.elapsed().user/nano; }
    double getSystemTime() const { return t.elapsed().system/nano; }*/
    void start() {
      clock_start= std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_seconds=std::chrono::system_clock::now()-clock_start;
      time_used= (double) elapsed_seconds.count();
    }
    void restart() { clock_start=std::chrono::system_clock::now();}
    void stop() {
      std::chrono::duration<double> elapsed_seconds=std::chrono::system_clock::now()-clock_start;
      time_used =  time_used + (double) elapsed_seconds.count(); }
    double getWallTime() const { return (double) time_used; }
    double getUserTime() const { return (double) time_used; }
    double getSystemTime() const { return (double) time_used; }
   //  double gettime(){
      //      return MPI_Wtime();
      // timespec time1;
      // clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
      // return (double) time1.tv_sec+1.0E-9*time1.tv_nsec;
    //   return (double) time_used;
    // }
    friend std::ostream& operator<<(std::ostream &o, const Timer &time_used) {
        int old_prec;
        GET_PRINT_PRECISION(old_prec);
        SET_PRINT_PRECISION(3);
        o << "    user   " << std::setw(6) << time_used;
        o << "    wall   " << std::setw(6) << time_used;
        SET_PRINT_PRECISION(old_prec);
        return o;
	}
private:
    //    boost::timer::cpu_timer t;
    double time_used;
    //std::chrono::duration<double> time_used;
    std::chrono::system_clock::time_point clock_start;
};
#endif // TIMER_H
