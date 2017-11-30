/*
 *  \date Jul 24, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 */

#include <fstream>
#include "Printer.h"
#include "Timer.h"
#include "parallel.h"

using namespace std;

static ofstream tp_outfile;
int Printer::printLevel = -1;
int Printer::precision = 12;
std::ostream *Printer::out = &std::cout;

void Printer::init(int level, const char *file) {
    int mpi_rank = mpiOrbRank;
    int mpi_size = mpiOrbSize;
    if (file != 0) {
        stringstream fname;
        if (mpi_size > 1) {
            fname << file << "-" << mpi_rank << ".out";
        } else {
            fname << file << ".out";
        }
        if (not fname.str().empty()) {
            tp_outfile.open(fname.str().c_str());
            setOutputStream(tp_outfile);
        }
    } else {
        if (mpi_rank > 0) {
            setPrintLevel(-1); // Higher ranks be quiet
        }
    }
    setScientific();
    setPrintLevel(level);
}

void Printer::printEnvironment(int level) {
    printout(level, endl);
    printSeparator(level, '-', 1);
    println(level, " MRCPP version   : " << PROGRAM_VERSION);
    println(level, " Git revision    : " << GIT_REVISION << endl);
    println(level, " Print level     : " << getPrintLevel());

#ifdef HAVE_BLAS
    println(level, " Linear algebra  : BLAS");
#else
    println(level, " Linear algebra  : EIGEN");
#endif

    int nHosts = mpiOrbSize;
    int nThreads = omp_get_max_threads();
#ifdef HAVE_MPI
#ifdef HAVE_OPENMP
    println(level, " Parallelization : MPI/OpenMP");
    println(level, " - MPI hosts     : " << nHosts);
    println(level, " - OMP threads   : " << nThreads);
    println(level, " - Total cores   : " << nThreads*nHosts);
#else
    println(level, " Parallelization : MPI");
    println(level, " - MPI hosts     : " << nHosts);
#endif
#else
#ifdef HAVE_OPENMP
    println(level, " Parallelization : OpenMP");
    println(level, " - OMP threads   : " << nThreads);
#else
    println(level, " Parallelization : NONE");
#endif
#endif

    printout(level, endl);
    printSeparator(level, '-', 2);
}


void Printer::printSeparator(int level, const char &sep, int newlines) {
    int N = 60;
    for (int i = 0; i < N; i++) {
        printout(level, sep);
    }
    for (int i = 0; i <= newlines; i++) {
        printout(level, endl);
    }
}

void Printer::printHeader(int level, const string &str, int newlines) {
    int N = 60;
    int len = str.size();
    printSeparator(level, '=', 0);
    int spaces = (N - len)/2;
    for (int i = 0; i < spaces; i++) {
        printout(level, " ");
    }
    println(level, str);
    printSeparator(level, '-', newlines);
}

void Printer::printFooter(int level, const Timer &t, int newlines) {
    printSeparator(level, '-');
    println(level, "                 Wall time: " << setw(11) << t << " sec ");
    printSeparator(level, '=', newlines);
}

void Printer::printDouble(int level, const std::string &name, double d) {
    char cName[31] = "                              ";
    for (int i = 0; i < 31; i++) {
        if (i < name.size()) {
            cName[i+1] = name[i];
        }
    }
    int oldPrec = Printer::setPrecision(5);
    println(level, cName << setw(29) << d);
    Printer::setPrecision(oldPrec);
}

void Printer::printTree(int level, const std::string &name, int n, double t) {
    char cName[31] = "                              ";
    for (int i = 0; i < 31; i++) {
        if (i < name.size()) {
            cName[i+1] = name[i];
        }
    }
    int oldPrec = Printer::setPrecision(5);
    println(level, cName << setw(12) << n << setw(17) << t);
    Printer::setPrecision(oldPrec);
}

int Printer::setPrintLevel(int i) {
    int oldLevel = printLevel;
    printLevel = i;
    return oldLevel;
}

int Printer::setPrecision(int i) {
    int oldPrec = precision;
    precision = i;
    *out << std::setprecision(i);
    return oldPrec;
}
