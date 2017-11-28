#pragma once

#include "config.h"

#define EIGEN_DONT_PARALLELIZE

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef HAVE_OPENMP
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

template<int D> class FunctionTree;
template<int D> class SerialFunctionTree;

namespace mpi {

#ifdef HAVE_MPI

template<int D>
struct tree_request {
    tree_request(FunctionTree<D> &tree, int t) {
        this->tag = t;
        this->s_tree = tree.getSerialFunctionTree();
        this->n_chunks = this->s_tree->nodeChunks.size();
        this->mpi_req = MPI_REQUEST_NULL;
    }
    int tag;
    int n_chunks;
    MPI_Status mpi_stat;
    MPI_Request mpi_req;
    SerialFunctionTree<D> *s_tree;
};

class communicator {
public:
    communicator() { }
    communicator &operator=(MPI_Comm c) { this->comm = c; return *this; }

    int size() const { int s; MPI_Comm_size(this->comm, &s); return s; }
    int rank() const { int r; MPI_Comm_rank(this->comm, &r); return r; }
    void barrier() const { MPI_Barrier(this->comm); }

    template<int D> void send_tree(tree_request<D> &request, int dst);
    template<int D> void recv_tree(tree_request<D> &request, int src);

private:
    MPI_Comm comm;
};

#else

template<int D>
struct tree_request {
    tree_request(FunctionTree<D> &tree, int t) { }
};

class communicator {
public:
    communicator() { }
    int rank() const { return 0; }
    int size() const { return 1; }
    void barrier() const { }

    template<int D> void send_tree(tree_request<D> &request, int dst) { }
    template<int D> void recv_tree(tree_request<D> &request, int src) { }
};

#endif

extern mpi::communicator world;

};
