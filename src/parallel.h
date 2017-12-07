#pragma once

#include <vector>
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

class requests {
public:
    requests() { }

    int size() const { return mpi_req.size(); }
    MPI_Request& operator[](int i) { return mpi_req[i]; }

    void clear() { mpi_req.clear(); }
    void push_back(MPI_Request req) { mpi_req.push_back(req); }

private:
    std::vector<MPI_Request> mpi_req;
};

class communicator {
public:
    communicator() { }
    communicator &operator=(MPI_Comm c) { this->comm = c; return *this; }

    int size() const { int s; MPI_Comm_size(this->comm, &s); return s; }
    int rank() const { int r; MPI_Comm_rank(this->comm, &r); return r; }
    void barrier() const { MPI_Barrier(this->comm); }

    template<int D> void send_tree(FunctionTree<D> &tree, int dst, int tag);
    template<int D> void recv_tree(FunctionTree<D> &tree, int src, int tag);

    template<int D> void isend_tree(FunctionTree<D> &tree, int dst, int tag, mpi::requests &req);
    void wait(mpi::requests &req);

private:
    MPI_Comm comm;
};

#else

typedef int requests;

class communicator {
public:
    communicator() { }
    int rank() const { return 0; }
    int size() const { return 1; }
    void barrier() const { }

    template<int D> void send_tree(FunctionTree<D> &tree, int dst, int tag) { }
    template<int D> void recv_tree(FunctionTree<D> &tree, int src, int tag) { }

    template<int D> void isend_tree(FunctionTree<D> &tree, int dst, int tag, requests &reqs) { }
    void wait(requests &reqs) { }
};

#endif

extern mpi::communicator world;

};
