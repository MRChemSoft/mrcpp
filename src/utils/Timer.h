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

#pragma once

#include <chrono>
#include <iostream>

namespace mrcpp {

using timeT = std::chrono::time_point<std::chrono::high_resolution_clock>;

class Timer final {
public:
    Timer(bool start_timer = true);
    Timer(const Timer &timer);
    Timer &operator=(const Timer &timer);

    void start();
    void resume();
    void stop();

    double getWallTime() const;

    friend std::ostream &operator<<(std::ostream &o, const Timer &timer) { return timer.print(o); }

private:
    bool running{false};
    double time_used{0.0};
    timeT clock_start;

    timeT now();
    double diffTime(timeT t2, timeT t1);

    std::ostream &print(std::ostream &o) const;
};

} // namespace mrcpp
