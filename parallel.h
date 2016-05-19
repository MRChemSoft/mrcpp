/*
 *
 *
 *  \date Sep 27, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif Convenience definitions for parallel builds
 */

#ifndef PARALLEL_H_
#define PARALLEL_H_

#include <boost/serialization/serialization.hpp>

#include "config.h"

int get_locale_index_range(int rank, int nWork, int &start, int &end);
int get_locale_index(int nWork, int idx);
bool locale_needs_sync(int nWork);

#ifdef HAVE_OPENMP

#define OPENMP
#include <omp.h>

#else

#define omp_get_max_threads() 1
#define omp_get_num_threads() 1
#define omp_get_thread_num() 0
#define omp_set_dynamic(n)
#define omp_set_lock(x)
#define omp_unset_lock(x)
#define omp_test_lock(x)

#endif

#ifdef HAVE_MPI

#include <boost/mpi.hpp>
namespace mpi = boost::mpi;
#define BOOST_MPI_HOMOGENEOUS

#else

namespace mpi {
struct communicator {
    int rank() const { return 0; }
    int size() const { return 1; }
    void barrier() const { }
};
struct environment {
    environment(int argc, char **argv) { }
};
typedef int request;
}

#endif

extern mpi::communicator world;
extern mpi::communicator node_group;

//template<int D>
//void broadcast_index_list(NodeIndexSet &idxSet) {
//#ifdef HAVE_MPI
//    if (node_group.size() > 1) {
//        static std::vector<std::vector<NodeIndex<D> > > commIdx;
//        commIdx.clear();

//        std::vector<NodeIndex<D> > tmpList;
//        typename std::set<const NodeIndex<D> *>::iterator it;
//        for (it = idxSet.begin(); it != idxSet.end(); it++) {
//            tmpList.push_back(**it);
//        }
//        mpi::all_gather(node_group, tmpList, commIdx);
//        idxSet.clear();
//        for (unsigned int i = 0; i < commIdx.size(); i++) {
//            for (unsigned int j = 0; j < commIdx[i].size(); j++) {
//                idxSet.insert(&commIdx[i][j]);
//            }
//        }
//    }
//#endif
//}

#endif /* PARALLEL_H_ */
