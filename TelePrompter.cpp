/*
 *  \date Jul 24, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 */

#include <fstream>
#include "TelePrompter.h"
#include "Timer.h"
#include "parallel.h"

using namespace std;

static ofstream tp_outfile;
int TelePrompter::printLevel = 0;
int TelePrompter::precision = 12;
std::ostream *TelePrompter::out = &std::cout;

void TelePrompter::printSeparator(const char &sep, int newlines) {
    int N = 60;
    for (int i = 0; i < N; i++) {
        cout << sep;
    }
    for (int i = 0; i <= newlines; i++) {
        cout << endl;
    }
}
void TelePrompter::printHeader(const string &str, int newlines) {
    int N = 60;
    int len = str.size();
    printSeparator('=', 0);
    int spaces = (N - len)/2;
    for (int i = 0; i < spaces; i++) {
        cout << " ";
    }
    cout << str << endl;
    printSeparator('-', newlines);
}

void TelePrompter::printFooter(const Timer &t, int newlines) {
    printSeparator('-');
    cout << " Elapsed time:     " << t << endl;
    printSeparator('=', newlines);
}

void TelePrompter::init(int printLevel, bool teletype, const char *fil) {
    SET_PRINT_PRECISION(15);
    cout << scientific << setprecision(14);

    mpi::communicator world;
    int rank=world.rank();
    SET_PRINT_LEVEL(printLevel);
    if (teletype and fil != 0) {
        stringstream fname;
        if (rank == 0) {
            if (printLevel < 0) {
                if (world.size() > 1) {
                    fname << fil << "-" << rank << ".out";
                } else {
                    fname << fil << ".out";
                }
            }
        } else {
            fname << fil << "-" << rank << ".out";
        }

        if (printLevel < 0) {
            SET_PRINT_LEVEL(-printLevel);
        }
        if (not fname.str().empty()) {
            tp_outfile.open(fname.str().c_str());
            SET_MESSAGE_STREAM(tp_outfile);
        }
    } else {
        if (rank > 0) {
            SET_PRINT_LEVEL(-10);
        }
    }
}
