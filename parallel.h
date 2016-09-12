#ifndef PARALLEL_H_
#define PARALLEL_H_

#include "config.h"
#include "SerialTree.h"

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

#include "SerialTree.h"

extern int MPI_rank;
extern int MPI_size;

void define_groups();
void MPI_Initializations();


#ifdef HAVE_MPI

template<int D>
void SendRcv_SerialTree(FunctionTree<D>* Tree, int source, int dest, int tag, MPI_Comm comm);
//void Share_memory(MPI_Comm ncomm, MPI_Comm ncomm_sh, int sh_size, double * d_ptr);

#else


#endif
#endif /* PARALLEL_H_ */
