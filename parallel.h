#ifndef PARALLEL_H_
#define PARALLEL_H_

#include "config.h"

template<int D> class FunctionTree;
class Orbital;

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

extern int MPI_Orb_rank;
extern int MPI_Orb_size;
extern int MPI_SH_rank;
extern int MPI_SH_size;
extern int MPI_SH_group_rank;
extern int MPI_SH_group_size;


void define_MPI_groups();
void MPI_Initializations();


#ifdef HAVE_MPI
#include <mpi.h>

extern MPI_Comm MPI_Comm_Orb;
extern MPI_Comm MPI_Comm_SH;
extern MPI_Comm MPI_Comm_SH_group;


/** Share memory within a compute node
 */
class SharedMemory {
public:    
    SharedMemory(int sh_size);
    ~SharedMemory();
    void allocShmem(int sh_size);
    double * sh_start_ptr; //start of shared block
    double * sh_max_ptr; //end of shared block
    double * sh_end_ptr; //end of used part
#ifdef HAVE_MPI
    MPI_Win sh_win; //MPI window object 
#endif
};

template<int D>
void Send_SerialTree(FunctionTree<D>* Tree, int Nchunks, int dest, int tag, MPI_Comm comm);
template<int D>
void IRcv_SerialTree(FunctionTree<D>* Tree, int Nchunks, int source, int tag, MPI_Comm comm);
template<int D>
void ISend_SerialTree(FunctionTree<D>* Tree, int Nchunks, int dest, int tag, MPI_Comm comm, MPI_Request& request);
template<int D>
void Rcv_SerialTree(FunctionTree<D>* Tree, int Nchunks, int source, int tag, MPI_Comm comm);
void Assign_NxN(int N, int* doi, int*doj, int* sendto, int* sendorb, int* rcvorb, int* MaxIter);
void Assign_NxN_sym(int N, int* doi, int*doj, int* sendto, int* sendorb, int* rcvorb, int* MaxIter);

void Share_memory(int sh_size, double * &d_ptr, MPI_Win & MPI_win);

#else


#endif
#endif /* PARALLEL_H_ */
