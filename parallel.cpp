#include "parallel.h"

mpi::communicator world;
mpi::communicator node_group;

int get_locale_index_range(int rank, int nWork, int &start, int &end) {
#ifdef HAVE_MPI
    int nLocales = node_group.size();
    int nPerLocale = (nWork - 1) / nLocales;
    int nLeft = nWork - nPerLocale * nLocales;
    if (rank < nLeft) {
        start = rank * (nPerLocale + 1);
        end = start + (nPerLocale + 1);
    } else {
        start = nLeft + rank * nPerLocale;
        end = start + nPerLocale;
    }
    return end - start;
#else
    start = 0;
    end = nWork;
    return nWork;
#endif
}

/** Return the host index having work unit "idx" */
int get_locale_index(int nWork, int idx) {
#ifdef HAVE_MPI
    int nLocales = node_group.size();
    int nPerLocale = (nWork - 1) / nLocales;
    int nLeft = nWork - nPerLocale * nLocales;
    int start, end;
    int rank;
    for (rank = 0; rank < nLocales; rank++) {
        if (rank < nLeft) {
            start = rank * (nPerLocale + 1);
            end = start + (nPerLocale + 1);
        } else {
            start = nLeft + rank * nPerLocale;
            end = start + nPerLocale;
        }
        if (idx >= start and idx < end) {
            break;
        }
    }
    return rank;
#else
    return 0;
#endif
}

/** If the locale has done less work than others, we need to ensure that
 they still can sync data from this locale. */
bool locale_needs_sync(int nWork) {
#ifdef HAVE_MPI
    int rank = node_group.rank();
    int nLocales = node_group.size();
    int nPerLocale = (nWork - 1) / nLocales;
    int nLeft = nWork - nPerLocale * nLocales;
    if (rank >= nLeft) {
        return true;
    }
#endif
    return false;
}
