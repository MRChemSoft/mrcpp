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

class Printer final {
public:
    static void init(int level = 0, int rank = 0, int size = 1, const char *file = nullptr);

    static void setOutputStream(std::ostream &o) { out = &o; }
    static void setScientific() { *out << std::scientific; }
    static void setFixed() { *out << std::fixed; }

    static int setWidth(int i);
    static int setPrecision(int i);
    static int setPrintLevel(int i);

    static int getWidth() { return printWidth; }
    static int getPrecision() { return printPrec; }
    static int getPrintLevel() { return printLevel; }

    static std::ostream *out;

private:
    static int printWidth;
    static int printLevel;
    static int printPrec;
    static int printRank;
    static int printSize;
};

namespace print {
void environment(int level);
void separator(int level, const char &sep, int newlines = 0);
void header(int level, const std::string &txt, int newlines = 0, const char &c = '=');
void footer(int level, const Timer &timer, int newlines = 0, const char &c = '=');
void memory(int level, const std::string &txt);
void value(int level, const std::string &txt, double v, const std::string &unit = "", int p = -1, bool sci = true);
void time(int level, const std::string &txt, const Timer &timer);
void tree(int level, const std::string &txt, int n, int m, double t);
template <int D> void tree(int level, const std::string &txt, const MWTree<D> &tree, const Timer &timer);
} // namespace print

// clang-format off
#define println(level, STR)                                                                                                             \
    { if (level <= mrcpp::Printer::getPrintLevel()) *mrcpp::Printer::out << STR << std::endl; }
#define printout(level, STR)                                                                                                            \
    { if (level <= mrcpp::Printer::getPrintLevel()) *mrcpp::Printer::out << STR; }

#define MSG_INFO(X)                                                                                                                     \
    {                                                                                                                                   \
        *mrcpp::Printer::out << "Info: " << __FILE__ << ": " << __func__ << "(), line " << __LINE__ << ": " << X << std::endl;          \
    }
#define MSG_WARN(X)                                                                                                                     \
    {                                                                                                                                   \
        *mrcpp::Printer::out << "Warning: " << __func__ << "(), line " << __LINE__ << ": " << X << std::endl;                           \
    }
#define MSG_ERROR(X)                                                                                                                    \
    {                                                                                                                                   \
        *mrcpp::Printer::out << "Error: " << __func__ << "(), line " << __LINE__ << ": " << X << std::endl;                             \
    }
#define MSG_ABORT(X)                                                                                                                    \
    {                                                                                                                                   \
        *mrcpp::Printer::out << "Error: " << __FILE__ << ": " << __func__ << "(), line " << __LINE__ << ": " << X << std::endl;         \
        abort();                                                                                                                        \
    }
#define INVALID_ARG_ABORT                                                                                                               \
    {                                                                                                                                   \
        *mrcpp::Printer::out << "Error, invalid argument passed: " << __func__ << "(), line " << __LINE__ << std::endl;                 \
        abort();                                                                                                                        \
    }
#define NOT_IMPLEMENTED_ABORT                                                                                                           \
    {                                                                                                                                   \
        *mrcpp::Printer::out << "Error: Not implemented, " << __FILE__ ", " << __func__ << "(), line " << __LINE__ << std::endl;        \
        abort();                                                                                                                        \
    }
#define NOT_REACHED_ABORT                                                                                                               \
    {                                                                                                                                   \
        *mrcpp::Printer::out << "Error, should not be reached: " << __func__ << "(), line " << __LINE__ << std::endl;                   \
        abort();                                                                                                                        \
    }
#define NEEDS_TESTING                                                                                                                   \
    {                                                                                                                                   \
        static bool __once = true;                                                                                                      \
        if (__once) {                                                                                                                   \
            __once = false;                                                                                                             \
            *mrcpp::Printer::out << "NEEDS TESTING: " << __FILE__ << ", " << __func__ << "(), line " << __LINE__ << std::endl;          \
        }                                                                                                                               \
    }
#define NEEDS_FIX(X)                                                                                                                    \
    {                                                                                                                                   \
        static bool __once = true;                                                                                                      \
        if (__once) {                                                                                                                   \
            __once = false;                                                                                                             \
            *mrcpp::Printer::out << "NEEDS FIX: " << __FILE__ << ", " << __func__ << "(), line " << __LINE__ << ": " << X << std::endl; \                                                                        \
        }                                                                                                                               \
    }
// clang-format on

} // namespace mrcpp
