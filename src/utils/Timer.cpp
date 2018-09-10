#include "Timer.h"
#include "Printer.h"

namespace mrcpp {

Timer::Timer(bool start_timer) : running(false), time_used(0.0) {
    if (start_timer) {
        start();
    }
}

Timer& Timer::operator=(const Timer &timer) {
    if (this != &timer) {
        this->time_used = timer.time_used;
        this->clock_start = timer.clock_start;
        this->running = timer.running;
    }
    return *this;
}

void Timer::start() {
    this->clock_start = now();
    this->time_used = 0.0;
    this->running = true;
}

void Timer::resume() {
    if (this->running) MSG_WARN("Timer already running");
    this->clock_start = now();
    this->running = true;
}

void Timer::stop() {
    if (not this->running) MSG_WARN("Timer not running");
    this->time_used += diffTime(now(), this->clock_start);
    this->running = false;
}

double Timer::getWallTime() const {
    if (this->running) MSG_WARN("Timer still running");
    return this->time_used;
}

timeT Timer::now(){
    return std::chrono::high_resolution_clock::now();
}

double Timer::diffTime(timeT t2, timeT t1) {
    std::chrono::duration<double> diff = t2-t1;
    return diff.count();
}

std::ostream& Timer::print(std::ostream &o) const {
    int old_prec = Printer::setPrecision(5);
    o << getWallTime();
    Printer::setPrecision(old_prec);
    return o;
}

} // namespace mrcpp
