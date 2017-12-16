/** \file Printer.h
 * Collection of assertions and a standard error/warn/info/debug
 * message interface.
 *
 */
#pragma once

#include <cassert>
#include <sstream>
#include <iostream>
#include <string>
#include <iomanip>

class Timer;

class Printer {
public:
    static void init(int level = 0, int rank = 0, int size = 1, const char *file = 0);
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

    static std::ostream *out;
private:
    static int printLevel;
    static int printPrec;
    static int printRank;
    static int printSize;
};

#define STR_DEBUG(S,X) {std::ostringstream _str;\
    _str << "Debug: " << __func__ << ",  line " << __LINE__ << \
    " in  " << __FILE__ << ": " << X << std::endl; S=_str.str();}
#define STR_INFO(S,X) {std::ostringstream _str;\
    _str << "Info: " << __func__ << ",  line " << __LINE__ << \
    " in  " << __FILE__ << ": " << X << std::endl; S=_str.str();}
#define STR_WARN(S,X) {std::ostringstream _str;\
    _str << "Warning: " << __func__ << ",  line " << __LINE__ << \
    " in  " << __FILE__ << ": " << X << std::endl; S=_str.str();}
#define STR_ERROR(S,X) {std::ostringstream _str;\
    _str << "Error: " << __func__ << ",  line " << __LINE__ << \
    " in  " << __FILE__ << ": " << X << std::endl; S=_str.str();}

#define println(level, STR) if (level <= Printer::getPrintLevel()) \
*Printer::out << STR << std::endl;
#define printout(level, STR) if (level <= Printer::getPrintLevel()) \
*Printer::out << STR;

#define MSG_DEBUG(X) { *Printer::out << "Debug: " << \
    __func__ << "(), line " << __LINE__ << "in " << \
    __FILE__ << ": " << X << std::endl;}
#define MSG_INFO(X) { *Printer::out << "Info: " << __FILE__ << ": " <<\
    __func__ << "(), line " << __LINE__ << ": " << X << std::endl;}
#define MSG_NOTE(X) { *Printer::out << "Note: " << __FILE__ << ": " <<\
    __func__ << "(), line " << __LINE__ << ": " << X << std::endl;}
#define MSG_WARN(X) { *Printer::out << "Warning: " \
    << __func__ << "(), line " << __LINE__ << ": " << X << std::endl;}
#define MSG_ERROR(X) { *Printer::out << "Error: " << \
    __func__ << "(), line " << __LINE__ << ": " << X << std::endl;}
#define MSG_FATAL(X) { *Printer::out << "Error: " << __FILE__ << ": " << \
    __func__ << "(), line " << __LINE__ << ": " << X << std::endl; abort();}

#define MSG_INVALID_ARG(X) { *Printer::out << \
    "Error, invalid argument passed: " << __func__ << "(), line " <<\
    __LINE__ << ": " << X << std::endl;}
#define INVALID_ARG_ABORT { *Printer::out << \
    "Error, invalid argument passed: " << __func__ << "(), line " <<\
    __LINE__ << std::endl; abort();}
#define NOT_REACHED_ABORT { *Printer::out << \
    "Error, should not be reached: " << __func__ << "(), line " <<\
    __LINE__ << std::endl; abort();}
#define INTERNAL_INCONSISTENCY { *Printer::out << \
    "Internal inconsistency! You have found a bug: " << __func__ << "(), line " <<\
    __LINE__ << std::endl; abort();}

#define NEEDS_TESTING {static bool __once = true;\
    if (__once) { __once = false;\
    *Printer::out << "NEEDS TESTING: " << __FILE__ << ", " << \
    __func__ << "(), line " << __LINE__ <<  std::endl;}}

#define ASSERT_FILE(A,B) {if (A == NULL) {\
    *Printer::out << "Error: " << __func__ << "(), line " << __LINE__ << \
    ": No such file, " << B << std::endl; abort();}}

#define NOT_IMPLEMENTED_ABORT {*Printer::out << \
    "Error: Not implemented, " << __FILE__ ", " << __func__ << "(), line " << __LINE__ << \
    std::endl; abort();}

#define NOTE(X) {static bool __once = true;\
    if (__once) { __once = false; *Printer::out << \
    "NOTE: " << __FILE__ << ", " << __func__ << "(), line " << __LINE__ << ": " << X << \
    std::endl; }}

#define NEEDS_FIX(X) {static bool __once = true;\
    if (__once) { __once = false; *Printer::out << \
    "NEEDS FIX: " << __FILE__ << ", " << __func__ << "(), line " << __LINE__ << ": " << X << \
    std::endl; }}

#define WRONG(X) {*Printer::out << \
    "WRONG: " << __FILE__ << ", " << __func__ << "(), line " << __LINE__ << ": " << X << \
    std::endl; abort();}

#define STR_DEBUG(S,X) {std::ostringstream _str;\
    _str << "Debug: " << __func__ << ",  line " << __LINE__ << \
    " in  " << __FILE__ << ": " << X << std::endl; S=_str.str();}
#define STR_INFO(S,X) {std::ostringstream _str;\
    _str << "Info: " << __func__ << ",  line " << __LINE__ << \
    " in  " << __FILE__ << ": " << X << std::endl; S=_str.str();}
#define STR_WARN(S,X) {std::ostringstream _str;\
    _str << "Warning: " << __func__ << ",  line " << __LINE__ << \
    " in  " << __FILE__ << ": " << X << std::endl; S=_str.str();}
#define STR_ERROR(S,X) {std::ostringstream _str;\
    _str << "Error: " << __func__ << ",  line " << __LINE__ << \
    " in  " << __FILE__ << ": " << X << std::endl; S=_str.str();}
