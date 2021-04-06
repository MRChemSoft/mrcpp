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

/** \file Printer.h
 * Collection of assertions and a standard error/warn/info/debug
 * message interface.
 *
 */
#pragma once

#include <cassert>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

namespace mrcpp {

class Timer;
template <int D> class MWTree;

/** @class Printer
 *
 * @brief Convenience class to handle printed output
 *
 * @details The ``Printer`` singleton class holds the current state of the print
 * environment. All ``mrcpp::print`` functions, as well as the ``println`` and
 * ``printout`` macros, take an integer print level as first argument. When the
 * global ``mrcpp::Printer`` is initialized with a given print level, only print
 * statements with a *lower* print level will be displayed. All internal printing
 * in MRCPP is at print level 10 or higher, so there is some flexibility left
 * (levels 0 through 9) for adjusting the print volume within the host program.
 *
 */

class Printer final {
public:
    static void init(int level = 0, int rank = 0, int size = 1, const char *file = nullptr);

    /** @brief Use scientific floating point notation, e.g. 1.0e-2 */
    static void setScientific() { *out << std::scientific; }

    /** @brief Use fixed floating point notation, e.g. 0.01 */
    static void setFixed() { *out << std::fixed; }

    /** @brief Set new line width for printed output
     *  @param[in] i: New width (number of characters)
     *  @returns Old width (number of characters)
     */
    static int setWidth(int i) {
        int oldWidth = printWidth;
        printWidth = i;
        return oldWidth;
    }

    /** @brief Set new precision for floating point output
     *  @param[in] i: New precision (digits after comma)
     *  @returns Old precision (digits after comma)
     */
    static int setPrecision(int i) {
        int oldPrec = printPrec;
        printPrec = i;
        *out << std::setprecision(i);
        return oldPrec;
    }

    /** @brief Set new print level
     *  @param[in] i: New print level
     *  @returns Old print level
     */
    static int setPrintLevel(int i) {
        int oldLevel = printLevel;
        printLevel = i;
        return oldLevel;
    }

    /** @returns Current line width (number of characters) */
    static int getWidth() { return printWidth; }

    /** @returns Current precision for floating point output (digits after comma) */
    static int getPrecision() { return printPrec; }

    /** @returns Current print level */
    static int getPrintLevel() { return printLevel; }

    static std::ostream *out;

private:
    static int printWidth;
    static int printLevel;
    static int printPrec;
    static int printRank;
    static int printSize;

    Printer() = delete; // No instances of this class

    static void setOutputStream(std::ostream &o) { out = &o; }
};

namespace print {
void environment(int level);
void separator(int level, const char &c, int newlines = 0);
void header(int level, const std::string &txt, int newlines = 0, const char &c = '=');
void footer(int level, const Timer &timer, int newlines = 0, const char &c = '=');
void memory(int level, const std::string &txt);
void value(int level, const std::string &txt, double v, const std::string &unit = "", int p = -1, bool sci = true);
void time(int level, const std::string &txt, const Timer &timer);
void tree(int level, const std::string &txt, int n, int m, double t);
template <int D> void tree(int level, const std::string &txt, const MWTree<D> &tree, const Timer &timer);
} // namespace print

// clang-format off

/** @brief Print text at the given print level, with newline */
#define println(level, STR)                                                                                                               \
    { if (level <= mrcpp::Printer::getPrintLevel()) *mrcpp::Printer::out << STR << std::endl; }

/** @brief Print text at the given print level, without newline */
#define printout(level, STR)                                                                                                              \
    { if (level <= mrcpp::Printer::getPrintLevel()) *mrcpp::Printer::out << STR; }

/** @brief Print info message */
#define MSG_INFO(STR)                                                                                                                     \
    {                                                                                                                                     \
        *mrcpp::Printer::out << "Info: " << __FILE__ << ": " << __func__ << "(), line " << __LINE__ << ": " << STR << std::endl;          \
    }

/** @brief Print warning message */
#define MSG_WARN(STR)                                                                                                                     \
    {                                                                                                                                     \
        *mrcpp::Printer::out << "Warning: " << __func__ << "(), line " << __LINE__ << ": " << STR << std::endl;                           \
    }

/** @brief Print error message, no abort*/
#define MSG_ERROR(STR)                                                                                                                    \
    {                                                                                                                                     \
        *mrcpp::Printer::out << "Error: " << __func__ << "(), line " << __LINE__ << ": " << STR << std::endl;                             \
    }

/** @brief Print error message and abort */
#define MSG_ABORT(STR)                                                                                                                    \
    {                                                                                                                                     \
        *mrcpp::Printer::out << "Error: " << __FILE__ << ": " << __func__ << "(), line " << __LINE__ << ": " << STR << std::endl;         \
        abort();                                                                                                                          \
    }

/** @brief You have passed an invalid argument to a function */
#define INVALID_ARG_ABORT                                                                                                                 \
    {                                                                                                                                     \
        *mrcpp::Printer::out << "Error, invalid argument passed: " << __func__ << "(), line " << __LINE__ << std::endl;                   \
        abort();                                                                                                                          \
    }

/** @brief You have reached a point in the code that is not yet implemented */
#define NOT_IMPLEMENTED_ABORT                                                                                                             \
    {                                                                                                                                     \
        *mrcpp::Printer::out << "Error: Not implemented, " << __FILE__ ", " << __func__ << "(), line " << __LINE__ << std::endl;          \
        abort();                                                                                                                          \
    }

/** @brief You have reached a point that should not be reached, bug or inconsistency */
#define NOT_REACHED_ABORT                                                                                                                 \
    {                                                                                                                                     \
        *mrcpp::Printer::out << "Error, should not be reached: " << __func__ << "(), line " << __LINE__ << std::endl;                     \
        abort();                                                                                                                          \
    }

/** @brief You have reached an experimental part of the code, results cannot be trusted */
#define NEEDS_TESTING                                                                                                                     \
    {                                                                                                                                     \
        static bool __once = true;                                                                                                        \
        if (__once) {                                                                                                                     \
            __once = false;                                                                                                               \
            *mrcpp::Printer::out << "NEEDS TESTING: " << __FILE__ << ", " << __func__ << "(), line " << __LINE__ << std::endl;            \
        }                                                                                                                                 \
    }

/** @brief You have hit a known bug that is yet to be fixed, results cannot be trusted */
#define NEEDS_FIX(STR)                                                                                                                    \
    {                                                                                                                                     \
        static bool __once = true;                                                                                                        \
        if (__once) {                                                                                                                     \
            __once = false;                                                                                                               \
            *mrcpp::Printer::out << "NEEDS FIX: " << __FILE__ << ", " << __func__ << "(), line " << __LINE__ << ": " << STR << std::endl; \                                                                        \
        }                                                                                                                                 \
    }
// clang-format on

} // namespace mrcpp
