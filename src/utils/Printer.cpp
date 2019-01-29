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

#include <fstream>

#include "config.h"

#include "Printer.h"
#include "Timer.h"

namespace mrcpp {

static std::ofstream tp_outfile;
int Printer::printLevel = -1;
int Printer::printPrec = 12;
int Printer::printRank = 0;
int Printer::printSize = 1;
std::ostream *Printer::out = &std::cout;

void Printer::init(int level, int rank, int size, const char *file) {
    printLevel = level;
    printRank = rank;
    printSize = size;
    if (file != nullptr) {
        std::stringstream fname;
        if (printSize > 1) {
            fname << file << "-" << printRank << ".out";
        } else {
            fname << file << ".out";
        }
        if (not fname.str().empty()) {
            tp_outfile.open(fname.str().c_str());
            setOutputStream(tp_outfile);
        }
    } else {
        if (printRank > 0) {
            setPrintLevel(-1); // Higher ranks be quiet
        }
    }
    setScientific();
    setPrecision(printPrec);
}

void Printer::printEnvironment(int level) {
    printout(level, std::endl);
    printSeparator(level, '-', 1);
    println(level, " MRCPP version   : " << PROGRAM_VERSION);
    println(level, " Git revision    : " << GIT_REVISION << std::endl);

#ifdef HAVE_BLAS
    println(level, " Linear algebra  : BLAS");
#else
    println(level, " Linear algebra  : EIGEN");
#endif

#ifdef HAVE_MPI
#ifdef HAVE_OPENMP
    println(level, " Parallelization : MPI/OpenMP");
#else
    println(level, " Parallelization : MPI");
#endif
#else
#ifdef HAVE_OPENMP
    println(level, " Parallelization : OpenMP");
#else
    println(level, " Parallelization : NONE");
#endif
#endif

    printout(level, std::endl);
    printSeparator(level, '-', 2);
}

void Printer::printSeparator(int level, const char &sep, int newlines) {
    int N = 60;
    for (int i = 0; i < N; i++) { printout(level, sep); }
    for (int i = 0; i <= newlines; i++) { printout(level, std::endl); }
}

void Printer::printHeader(int level, const std::string &str, int newlines) {
    int N = 60;
    int len = str.size();
    printSeparator(level, '=', 0);
    int spaces = (N - len) / 2;
    for (int i = 0; i < spaces; i++) { printout(level, " "); }
    println(level, str);
    printSeparator(level, '-', newlines);
}

void Printer::printFooter(int level, const Timer &t, int newlines) {
    printSeparator(level, '-');
    println(level, "                 Wall time: " << std::setw(11) << t << " sec ");
    printSeparator(level, '=', newlines);
}

void Printer::printDouble(int level, const std::string &str, double d, int p) {
    char cStr[31] = "                              ";
    for (int i = 0; i < 31; i++) {
        if (i < str.size()) { cStr[i + 1] = str[i]; }
    }
    int oldPrec = getPrecision();
    if (p > 0) setPrecision(p);
    println(level, cStr << std::setw(29) << d);
    setPrecision(oldPrec);
}

void Printer::printTree(int level, const std::string &str, int n, double t) {
    char cStr[31] = "                              ";
    for (int i = 0; i < 31; i++) {
        if (i < str.size()) { cStr[i + 1] = str[i]; }
    }
    int oldPrec = setPrecision(5);
    println(level, cStr << std::setw(12) << n << std::setw(17) << t);
    setPrecision(oldPrec);
}

void Printer::printTime(int level, const std::string &str, const Timer &t) {
    char cStr[31] = "                              ";
    for (int i = 0; i < 31; i++) {
        if (i < str.size()) { cStr[i + 1] = str[i]; }
    }
    int oldPrec = setPrecision(5);
    println(level, cStr << std::setw(29) << t);
    setPrecision(oldPrec);
}

int Printer::setPrintLevel(int i) {
    int oldLevel = printLevel;
    printLevel = i;
    return oldLevel;
}

int Printer::setPrecision(int i) {
    int oldPrec = printPrec;
    printPrec = i;
    *out << std::setprecision(i);
    return oldPrec;
}

// parse a string and returns the nth integer number
int Printer::getVal(char *line, int n) {
    char *p = line;
    int len = 0;
    for (int i = 0; i < n - 1; i++) {
        // jump over n-1 first numbers
        while (*p < '0' || *p > '9') p++;
        while (*p >= '0' && *p <= '9') p++;
    }
    while (*p < '0' || *p > '9') p++;
    char *s = p;
    while (*s >= '0' && *s <= '9') {
        s++;
        line[len] = p[len];
        len++;
    }
    p[len] = 0;
    return atoi(p);
}

/** Prints (and returns) the current memory usage of this process
 */
int Printer::printMem(char *txt, bool silent) {
    FILE *file = fopen("/proc/self/statm", "r");
    int val = -1;
    char line[80];
    while (fgets(line, 80, file) != nullptr) {
        val = getVal(line, 6); // sixth number is data+stack in pages (4kB)
        if (not silent) std::cout << &txt << val * 4.0 / (1024.0 * 1024) << "GB" << std::endl;
    }
    fclose(file);
    return val;
}

} // namespace mrcpp
