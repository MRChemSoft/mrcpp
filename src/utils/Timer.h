/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2021 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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

namespace mrcpp {

/**
 * @file
 * @brief Wall-clock timing utilities.
 */

/**
 * @typedef timeT
 * @brief Timestamp type used by Timer.
 *
 * @details
 * Alias for `std::chrono::time_point<std::chrono::high_resolution_clock>`.
 * The actual clock may map to a platform-specific high-resolution source.
 * It typically offers sub-microsecond resolution but is not guaranteed to be
 * steady on all standard libraries (i.e., it may jump if the underlying clock
 * is adjusted). The @ref Timer class measures *wall* time, not CPU time.
 */
using timeT = std::chrono::time_point<std::chrono::high_resolution_clock>;

/**
 * @class Timer
 * @brief Lightweight wall-time stopwatch with start/resume/stop semantics.
 *
 * @details
 * `Timer` accumulates elapsed *wall* time across one or more running intervals.
 * A newly constructed timer can optionally start immediately.
 *
 * ### State machine
 * - **Stopped** (default when constructed with `start_timer == false`):
 *   - `elapsed()` returns the accumulated time.
 *   - `resume()` starts a new interval without clearing accumulated time.
 *   - `start()` clears accumulated time and starts fresh from zero.
 * - **Running** (default when constructed with `start_timer == true`):
 *   - `elapsed()` returns the live time since the most recent start/resume,
 *     ignoring previously accumulated time until `stop()` is called.
 *   - `stop()` ends the current interval and adds it to the accumulation.
 *
 * ### Characteristics
 * - Measures wall time (affected by system sleep/suspend).
 * - Very low overhead; suitable for inner-loop timing in most cases.
 * - Not thread-safe: do not share a single instance across threads without
 *   external synchronization.
 *
 * ### Example
 * @code{.cpp}
 * mrcpp::Timer t;                // starts immediately by default
 * // ... code section A ...
 * t.stop();
 * double a = t.elapsed();         // seconds for section A
 *
 * t.resume();
 * // ... code section B ...
 * t.stop();
 * double total = t.elapsed();     // seconds for A + B
 *
 * t.start();                      // reset and start from zero
 * // ... code section C ...
 * double live = t.elapsed();      // live time while running
 * @endcode
 */
class Timer final {
public:
    /**
     * @brief Construct a timer.
     * @param start_timer If true, the timer is started immediately with
     *                    accumulated time cleared to zero.
     *
     * @note Default is `true` for convenience; pass `false` to construct
     *       a stopped timer and control the first start explicitly.
     */
    Timer(bool start_timer = true);

    /**
     * @brief Copy constructor.
     * @details Copies the running state, accumulated time, and the last start
     *          timestamp. If the source is running, the copy will also be
     *          running and will measure from the same start instant.
     */
    Timer(const Timer &timer);

    /**
     * @brief Copy assignment.
     * @details Assigns running state, accumulated time, and start timestamp.
     *          Self-assignment is a no-op.
     * @return Reference to `*this`.
     */
    Timer &operator=(const Timer &timer);

    /**
     * @brief Start from zero.
     * @details Resets the accumulated time to 0 and begins a new running
     *          interval starting "now". Use this to time a fresh region.
     */
    void start();

    /**
     * @brief Resume without clearing accumulated time.
     * @details If the timer is stopped, begins a new running interval starting
     *          "now". If already running, the call has no effect besides
     *          potentially issuing a diagnostic in the implementation.
     */
    void resume();

    /**
     * @brief Stop and accumulate.
     * @details Ends the current running interval and adds its duration to the
     *          accumulated time. If already stopped, the call has no effect
     *          besides potentially issuing a diagnostic in the implementation.
     */
    void stop();

    /**
     * @brief Get elapsed time in seconds.
     * @details
     * - If the timer is **running**, returns the time since the most recent
     *   `start()` or `resume()` (not including previously accumulated time).
     * - If the timer is **stopped**, returns the total accumulated time across
     *   all completed intervals since the last `start()`.
     *
     * @return Elapsed wall time in seconds as a `double`.
     */
    double elapsed() const;

private:
    /// @brief True if the timer is currently running.
    bool running{false};

    /// @brief Accumulated time in seconds across completed intervals.
    double time_used{0.0};

    /// @brief Timestamp when the current interval was started/resumed.
    timeT clock_start;

    /**
     * @brief Current timestamp helper.
     * @return `high_resolution_clock::now()`.
     */
    timeT now() const;

    /**
     * @brief Difference between two timestamps.
     * @param t2 Later timestamp.
     * @param t1 Earlier timestamp.
     * @return `(t2 - t1)` expressed in seconds as a `double`.
     */
    double diffTime(timeT t2, timeT t1) const;
};

} // namespace mrcpp