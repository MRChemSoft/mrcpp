/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2019 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
 *
 * This file is part of MRCPP.
 *
 * MRCPP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRCPP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRCPP.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRCPP, see:
 * <https://mrcpp.readthedocs.io/>
 */

#include "Timer.h"
#include "Printer.h"

namespace mrcpp {

Timer::Timer(bool start_timer) {
    if (start_timer) { start(); }
}

Timer::Timer(const Timer &timer)
        : running(timer.running)
        , time_used(timer.time_used)
        , clock_start(timer.clock_start) {}

Timer &Timer::operator=(const Timer &timer) {
    if (this != &timer) {
        this->running = timer.running;
        this->time_used = timer.time_used;
        this->clock_start = timer.clock_start;
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

timeT Timer::now() {
    return std::chrono::high_resolution_clock::now();
}

double Timer::diffTime(timeT t2, timeT t1) {
    std::chrono::duration<double> diff = t2 - t1;
    return diff.count();
}

std::ostream &Timer::print(std::ostream &o) const {
    int old_prec = Printer::setPrecision(5);
    o << getWallTime();
    Printer::setPrecision(old_prec);
    return o;
}

} // namespace mrcpp
