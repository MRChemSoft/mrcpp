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

class Printer final {
public:
    static void init(int level = 0, int rank = 0, int size = 1, const char *file = nullptr);
    static void printEnvironment(int level = 0);

    static void printSeparator(int level, const char &sep, int newlines = 0);
    static void printHeader(int level, const std::string &str, int newlines = 0);
    static void printFooter(int level, const Timer &t, int newlines = 0);

    static void printDouble(int level, const std::string &str, double d, int p = -1);
    static void printTree(int level, const std::string &str, int n, double t);
    static void printTime(int level, const std::string &str, const Timer &t);

    static void setOutputStream(std::ostream &o) { out = &o; }
    static void setScientific() { *out << std::scientific; }
    static void setFixed() { *out << std::fixed; }

    static int setPrecision(int i);
    static int setPrintLevel(int i);
    static int getPrecision() { return printPrec; }
    static int getPrintLevel() { return printLevel; }
    static int getVal(char *line, int n = 1);
    static int printMem(char *txt, bool silent = false);

    static std::ostream *out;

private:
    static int printLevel;
    static int printPrec;
    static int printRank;
    static int printSize;
};

#define STR_DEBUG(S, X)                                                                                                \
    {                                                                                                                  \
        std::ostringstream _str;                                                                                       \
        _str << "Debug: " << __func__ << ",  line " << __LINE__ << " in  " << __FILE__ << ": " << X << std::endl;      \
        S = _str.str();                                                                                                \
    }
#define STR_INFO(S, X)                                                                                                 \
    {                                                                                                                  \
        std::ostringstream _str;                                                                                       \
        _str << "Info: " << __func__ << ",  line " << __LINE__ << " in  " << __FILE__ << ": " << X << std::endl;       \
        S = _str.str();                                                                                                \
    }
#define STR_WARN(S, X)                                                                                                 \
    {                                                                                                                  \
        std::ostringstream _str;                                                                                       \
        _str << "Warning: " << __func__ << ",  line " << __LINE__ << " in  " << __FILE__ << ": " << X << std::endl;    \
        S = _str.str();                                                                                                \
    }
#define STR_ERROR(S, X)                                                                                                \
    {                                                                                                                  \
        std::ostringstream _str;                                                                                       \
        _str << "Error: " << __func__ << ",  line " << __LINE__ << " in  " << __FILE__ << ": " << X << std::endl;      \
        S = _str.str();                                                                                                \
    }

#define println(level, STR)                                                                                            \
    if (level <= mrcpp::Printer::getPrintLevel()) *mrcpp::Printer::out << STR << std::endl;
#define printout(level, STR)                                                                                           \
    if (level <= mrcpp::Printer::getPrintLevel()) *mrcpp::Printer::out << STR;

#define MSG_DEBUG(X)                                                                                                   \
    {                                                                                                                  \
        *mrcpp::Printer::out << "Debug: " << __func__ << "(), line " << __LINE__ << "in " << __FILE__ << ": " << X     \
                             << std::endl;                                                                             \
    }
#define MSG_INFO(X)                                                                                                    \
    {                                                                                                                  \
        *mrcpp::Printer::out << "Info: " << __FILE__ << ": " << __func__ << "(), line " << __LINE__ << ": " << X       \
                             << std::endl;                                                                             \
    }
#define MSG_NOTE(X)                                                                                                    \
    {                                                                                                                  \
        *mrcpp::Printer::out << "Note: " << __FILE__ << ": " << __func__ << "(), line " << __LINE__ << ": " << X       \
                             << std::endl;                                                                             \
    }
#define MSG_WARN(X)                                                                                                    \
    { *mrcpp::Printer::out << "Warning: " << __func__ << "(), line " << __LINE__ << ": " << X << std::endl; }
#define MSG_ERROR(X)                                                                                                   \
    { *mrcpp::Printer::out << "Error: " << __func__ << "(), line " << __LINE__ << ": " << X << std::endl; }
#define MSG_FATAL(X)                                                                                                   \
    {                                                                                                                  \
        *mrcpp::Printer::out << "Error: " << __FILE__ << ": " << __func__ << "(), line " << __LINE__ << ": " << X      \
                             << std::endl;                                                                             \
        abort();                                                                                                       \
    }

#define MSG_INVALID_ARG(X)                                                                                             \
    {                                                                                                                  \
        *mrcpp::Printer::out << "Error, invalid argument passed: " << __func__ << "(), line " << __LINE__ << ": " << X \
                             << std::endl;                                                                             \
    }
#define INVALID_ARG_ABORT                                                                                              \
    {                                                                                                                  \
        *mrcpp::Printer::out << "Error, invalid argument passed: " << __func__ << "(), line " << __LINE__              \
                             << std::endl;                                                                             \
        abort();                                                                                                       \
    }
#define NOT_REACHED_ABORT                                                                                              \
    {                                                                                                                  \
        *mrcpp::Printer::out << "Error, should not be reached: " << __func__ << "(), line " << __LINE__ << std::endl;  \
        abort();                                                                                                       \
    }
#define INTERNAL_INCONSISTENCY                                                                                         \
    {                                                                                                                  \
        *mrcpp::Printer::out << "Internal inconsistency! You have found a bug: " << __func__ << "(), line "            \
                             << __LINE__ << std::endl;                                                                 \
        abort();                                                                                                       \
    }

#define NEEDS_TESTING                                                                                                  \
    {                                                                                                                  \
        static bool __once = true;                                                                                     \
        if (__once) {                                                                                                  \
            __once = false;                                                                                            \
            *mrcpp::Printer::out << "NEEDS TESTING: " << __FILE__ << ", " << __func__ << "(), line " << __LINE__       \
                                 << std::endl;                                                                         \
        }                                                                                                              \
    }

#define ASSERT_FILE(A, B)                                                                                              \
    {                                                                                                                  \
        if (A == NULL) {                                                                                               \
            *mrcpp::Printer::out << "Error: " << __func__ << "(), line " << __LINE__ << ": No such file, " << B        \
                                 << std::endl;                                                                         \
            abort();                                                                                                   \
        }                                                                                                              \
    }

#define NOT_IMPLEMENTED_ABORT                                                                                          \
    {                                                                                                                  \
        *mrcpp::Printer::out << "Error: Not implemented, " << __FILE__ ", " << __func__ << "(), line " << __LINE__     \
                             << std::endl;                                                                             \
        abort();                                                                                                       \
    }

#define NOTE(X)                                                                                                        \
    {                                                                                                                  \
        static bool __once = true;                                                                                     \
        if (__once) {                                                                                                  \
            __once = false;                                                                                            \
            *mrcpp::Printer::out << "NOTE: " << __FILE__ << ", " << __func__ << "(), line " << __LINE__ << ": " << X   \
                                 << std::endl;                                                                         \
        }                                                                                                              \
    }

#define NEEDS_FIX(X)                                                                                                   \
    {                                                                                                                  \
        static bool __once = true;                                                                                     \
        if (__once) {                                                                                                  \
            __once = false;                                                                                            \
            *mrcpp::Printer::out << "NEEDS FIX: " << __FILE__ << ", " << __func__ << "(), line " << __LINE__ << ": "   \
                                 << X << std::endl;                                                                    \
        }                                                                                                              \
    }

#define WRONG(X)                                                                                                       \
    {                                                                                                                  \
        *mrcpp::Printer::out << "WRONG: " << __FILE__ << ", " << __func__ << "(), line " << __LINE__ << ": " << X      \
                             << std::endl;                                                                             \
        abort();                                                                                                       \
    }

#define STR_DEBUG(S, X)                                                                                                \
    {                                                                                                                  \
        std::ostringstream _str;                                                                                       \
        _str << "Debug: " << __func__ << ",  line " << __LINE__ << " in  " << __FILE__ << ": " << X << std::endl;      \
        S = _str.str();                                                                                                \
    }
#define STR_INFO(S, X)                                                                                                 \
    {                                                                                                                  \
        std::ostringstream _str;                                                                                       \
        _str << "Info: " << __func__ << ",  line " << __LINE__ << " in  " << __FILE__ << ": " << X << std::endl;       \
        S = _str.str();                                                                                                \
    }
#define STR_WARN(S, X)                                                                                                 \
    {                                                                                                                  \
        std::ostringstream _str;                                                                                       \
        _str << "Warning: " << __func__ << ",  line " << __LINE__ << " in  " << __FILE__ << ": " << X << std::endl;    \
        S = _str.str();                                                                                                \
    }
#define STR_ERROR(S, X)                                                                                                \
    {                                                                                                                  \
        std::ostringstream _str;                                                                                       \
        _str << "Error: " << __func__ << ",  line " << __LINE__ << " in  " << __FILE__ << ": " << X << std::endl;      \
        S = _str.str();                                                                                                \
    }

} // namespace mrcpp
