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
/**
 * @file
 * @brief Lightweight, level-based printing, diagnostics, and assertion helpers.
 *
 * This header provides:
 *  - A singleton-style @ref mrcpp::Printer that controls global print level,
 *    numeric formatting (scientific vs. fixed), width and precision, and the
 *    active output stream (stdout or a per-rank file).
 *  - A small set of convenience functions in @ref mrcpp::print for consistent,
 *    nicely formatted environment, headers/footers, timers, memory usage, and
 *    tree statistics.
 *  - A family of macros (e.g. @ref println, @ref MSG_ABORT) to emit messages
 *    with source context and optional termination semantics.
 *
 * ### Print-level convention
 * Every printing API takes an integer *level*. Output is produced iff
 * `level <= Printer::getPrintLevel()`. Internal MRCPP prints use \>= 10,
 * leaving levels 0–9 available for host/user code control.
 *
 * @warning The facilities herein are process-local; when used in MPI programs,
 *          each rank will emit messages independently unless explicitly gated
 *          by rank logic in the caller.
 */

#include <cassert>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

namespace mrcpp {

class Timer;
template <int D, typename T> class MWTree;

/**
 * @class Printer
 * @brief Process-local controller for formatted, level-gated output.
 *
 * @details
 * The Printer class is used in a singleton-like fashion via static methods.
 * Call init once near program start (optionally per MPI rank) to set:
 *  - global print level (messages with a higher level are suppressed),
 *  - rank/size metadata (used to route output),
 *  - destination stream (stdout or a per-rank file).
 *
 * After initialization, helper functions/macros (see mrcpp::print and the
 * macros below) consult the configured level and stream.
 *
 * @par Example
 * @code{.cpp}
 * // Rank 0 prints to screen, others remain silent:
 * Printer::init(5, rank, size);
 * println(2, "Hello at level 2");    // shown when level >= 2
 * println(8, "Debug details...");    // suppressed when level < 8
 * @endcode
 */
class Printer final {
public:
    /**
     * @brief Initialize printing environment.
     *
     * @param level   Maximum verbosity to emit (inclusive).
     * @param rank    MPI rank of this process (default 0).
     * @param size    MPI world size (default 1).
     * @param file    Optional base filename. If provided and @p size>1,
     *                output is written to "<file>-<rank>.out"; otherwise to
     *                "<file>.out". If null, output goes to stdout. When
     *                @p file is null and @p rank>0, printing is disabled
     *                by setting the print level to -1.
     *
     * @note Also sets scientific notation as the default numeric format.
     */
    static void init(int level = 0, int rank = 0, int size = 1, const char *file = nullptr);

    /** @brief Use scientific floating-point notation (e.g., 1.23e-2). */
    static void setScientific() { *out << std::scientific; }

    /** @brief Use fixed floating-point notation (e.g., 0.0123). */
    static void setFixed() { *out << std::fixed; }

    /**
     * @brief Set line width for formatted helpers.
     * @param i New width in characters.
     * @return Previous width.
     */
    static int setWidth(int i) {
        int oldWidth = printWidth;
        printWidth = i;
        return oldWidth;
    }

    /**
     * @brief Set precision for floating-point output.
     * @param i Digits after the decimal point.
     * @return Previous precision.
     */
    static int setPrecision(int i) {
        int oldPrec = printPrec;
        printPrec = i;
        *out << std::setprecision(i);
        return oldPrec;
    }

    /**
     * @brief Set the global print level threshold.
     * @param i New level; only messages with level <= @p i are printed.
     * @return Previous print level.
     */
    static int setPrintLevel(int i) {
        int oldLevel = printLevel;
        printLevel = i;
        return oldLevel;
    }

    /** @return Current line width (characters). */
    static int getWidth() { return printWidth; }

    /** @return Current floating-point precision (digits after decimal). */
    static int getPrecision() { return printPrec; }

    /** @return Current global print level threshold. */
    static int getPrintLevel() { return printLevel; }

    /**
     * @brief Active output stream (stdout or a file).
     * @warning Pointer is owned externally; @ref init manages it appropriately.
     */
    static std::ostream *out;

private:
    static int printWidth; ///< Line width used by @ref mrcpp::print helpers.
    static int printLevel; ///< Global verbosity threshold.
    static int printPrec;  ///< Floating-point precision (digits).
    static int printRank;  ///< MPI rank (for routing decisions).
    static int printSize;  ///< MPI world size.

    Printer() = delete; ///< Non-instantiable utility.

    /// @brief Redirect all output to @p o (used internally by @ref init).
    static void setOutputStream(std::ostream &o) { out = &o; }
};

/**
 * @namespace mrcpp::print
 * @brief Nicely formatted, level-aware printing helpers.
 *
 * These helpers produce standardized, aligned, and labeled output for common
 * diagnostics: environment summaries, section headers/footers, timers, memory
 * usage, and tree statistics. Each function is level-gated via
 * @ref Printer::getPrintLevel().
 */
namespace print {

/**
 * @brief Print MRCPP and build environment information.
 * @param level Activation level threshold.
 *
 * @details Includes library version, Git metadata, linear algebra backend,
 * and parallelization mode (MPI/OpenMP).
 */
void environment(int level);

/**
 * @brief Print a full separator line composed of @p c characters.
 * @param level Activation level.
 * @param c     Filler character (e.g., '-', '=').
 * @param newlines Number of extra trailing blank lines (default 0).
 */
void separator(int level, const char &c, int newlines = 0);

/**
 * @brief Print a centered header with a framed title.
 * @param level Activation level.
 * @param txt   Header text.
 * @param newlines Extra trailing blank lines (default 0).
 * @param c     Filler character for the frame (default '=').
 */
void header(int level, const std::string &txt, int newlines = 0, const char &c = '=');

/**
 * @brief Print a footer containing elapsed wall time and a closing frame.
 * @param level Activation level.
 * @param timer Timer whose @ref Timer::elapsed is shown.
 * @param newlines Extra trailing blank lines (default 0).
 * @param c     Filler character for the closing line (default '=').
 */
void footer(int level, const Timer &timer, int newlines = 0, const char &c = '=');

/**
 * @brief Print current process memory usage.
 * @param level Activation level.
 * @param txt   Label to show before the value (aligned).
 */
void memory(int level, const std::string &txt);

/**
 * @brief Print a labeled scalar with unit in aligned columns.
 * @param level Activation level.
 * @param txt   Label.
 * @param v     Value.
 * @param unit  Unit string (optional).
 * @param p     Precision; if negative uses @ref Printer::getPrecision.
 * @param sci   Scientific formatting when true, fixed when false.
 */
void value(int level, const std::string &txt, double v, const std::string &unit = "", int p = -1, bool sci = true);

/**
 * @brief Print an elapsed time value from a @ref Timer.
 * @param level Activation level.
 * @param txt   Label.
 * @param timer Timer whose @ref Timer::elapsed is shown.
 */
void time(int level, const std::string &txt, const Timer &timer);

/**
 * @brief Print tree statistics (nodes, memory, wall time).
 * @param level Activation level.
 * @param txt   Label/section title.
 * @param n     Number of nodes.
 * @param m     Memory usage in kB.
 * @param t     Elapsed wall time in seconds.
 */
void tree(int level, const std::string &txt, int n, int m, double t);

/**
 * @brief Print tree statistics extracted from an @ref MWTree and a @ref Timer.
 * @tparam D Spatial dimension.
 * @tparam T Coefficient scalar type.
 * @param level Activation level.
 * @param txt   Label/section title.
 * @param tree  Tree whose node/memory info is reported.
 * @param timer Timer whose elapsed time is reported.
 */
template <int D, typename T>
void tree(int level, const std::string &txt, const MWTree<D, T> &tree, const Timer &timer);

} // namespace print

// ============================================================================
// Macros : level-aware printing and diagnostics
// ============================================================================

/**
 * @def println(level, STR)
 * @brief Print text followed by newline if @p level is enabled.
 */
#define println(level, STR)                                                                       \
    { if (level <= mrcpp::Printer::getPrintLevel()) *mrcpp::Printer::out << STR << std::endl; }

/**
 * @def printout(level, STR)
 * @brief Print text without newline if @p level is enabled.
 */
#define printout(level, STR)                                                                      \
    { if (level <= mrcpp::Printer::getPrintLevel()) *mrcpp::Printer::out << STR; }

/**
 * @def MSG_INFO(STR)
 * @brief Emit an informational message with source location.
 */
#define MSG_INFO(STR)                                                                             \
    {                                                                                             \
        *mrcpp::Printer::out << "Info: " << __FILE__ << ": " << __func__                          \
                             << "(), line " << __LINE__ << ": " << STR << std::endl;              \
    }

/**
 * @def MSG_WARN(STR)
 * @brief Emit a warning message with source location.
 */
#define MSG_WARN(STR)                                                                             \
    {                                                                                             \
        *mrcpp::Printer::out << "Warning: " << __func__ << "(), line " << __LINE__                \
                             << ": " << STR << std::endl;                                         \
    }

/**
 * @def MSG_ERROR(STR)
 * @brief Emit a non-fatal error message with source location.
 */
#define MSG_ERROR(STR)                                                                            \
    {                                                                                             \
        *mrcpp::Printer::out << "Error: " << __func__ << "(), line " << __LINE__                  \
                             << ": " << STR << std::endl;                                         \
    }

/**
 * @def MSG_ABORT(STR)
 * @brief Emit an error message with source location and abort the process.
 */
#define MSG_ABORT(STR)                                                                            \
    {                                                                                             \
        *mrcpp::Printer::out << "Error: " << __FILE__ << ": " << __func__                         \
                             << "(), line " << __LINE__ << ": " << STR << std::endl;              \
        abort();                                                                                  \
    }

/**
 * @def INVALID_ARG_ABORT
 * @brief Abort with a standardized message for invalid arguments.
 */
#define INVALID_ARG_ABORT                                                                         \
    {                                                                                             \
        *mrcpp::Printer::out << "Error, invalid argument passed: " << __func__                    \
                             << "(), line " << __LINE__ << std::endl;                             \
        abort();                                                                                  \
    }

/**
 * @def NOT_IMPLEMENTED_ABORT
 * @brief Abort with a standardized message for unimplemented code paths.
 */
#define NOT_IMPLEMENTED_ABORT                                                                     \
    {                                                                                             \
        *mrcpp::Printer::out << "Error: Not implemented, " << __FILE__ << ", " << __func__        \
                             << "(), line " << __LINE__ << std::endl;                             \
        abort();                                                                                  \
    }

/**
 * @def NOT_REACHED_ABORT
 * @brief Abort for code paths that should be logically unreachable.
 */
#define NOT_REACHED_ABORT                                                                         \
    {                                                                                             \
        *mrcpp::Printer::out << "Error, should not be reached: " << __func__                      \
                             << "(), line " << __LINE__ << std::endl;                             \
        abort();                                                                                  \
    }

/**
 * @def NEEDS_TESTING
 * @brief Emit a one-time notice that a code path is experimental.
 *
 * Prints exactly once per process at the first hit, then stays quiet.
 */
#define NEEDS_TESTING                                                                             \
    {                                                                                             \
        static bool __once = true;                                                                \
        if (__once) {                                                                             \
            __once = false;                                                                       \
            *mrcpp::Printer::out << "NEEDS TESTING: " << __FILE__ << ", " << __func__             \
                                 << "(), line " << __LINE__ << std::endl;                         \
        }                                                                                         \
    }

/**
 * @def NEEDS_FIX(STR)
 * @brief Emit a one-time notice that a known issue affects this code path.
 * @param STR Short description of the known issue.
 */
#define NEEDS_FIX(STR)                                                                            \
    {                                                                                             \
        static bool __once = true;                                                                \
        if (__once) {                                                                             \
            __once = false;                                                                       \
            *mrcpp::Printer::out << "NEEDS FIX: " << __FILE__ << ", " << __func__                 \
                                 << "(), line " << __LINE__ << ": " << STR << std::endl;          \
        }                                                                                         \
    }

} // namespace mrcpp