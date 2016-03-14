/*
 *  \date Jul 24, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 */

#include <fstream>
#include "TelePrompter.h"
#include "parallel.h"

using namespace std;

static ofstream tp_outfile;
int TelePrompter::printLevel = 0;
int TelePrompter::precision = 12;
std::ostream *TelePrompter::out = &std::cout;

#define PRINT_SEPARATOR(A,N) for (int i=0; i < N; i++) cout << A; cout << endl;

void TelePrompter::printHeader(const std::string &str) {
    int len = 78;
    if (str.size() == 0 and str.size() < 74) {
        len = str.size() + 5;
    }
    PRINT_SEPARATOR("=", len);
    if (str.size() != 0) {
        cout << str << endl;
        PRINT_SEPARATOR("=", len);
    }
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
