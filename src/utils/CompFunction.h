#pragma once

#include "trees/FunctionTree.h"

using namespace Eigen;

namespace mrcpp {

template <int D>
struct CompFunctionData {
    // additional data that describe the overall multicomponent function (defined by user):
    // occupancy, quantum number, norm, etc.
    int Ncomp{1}; // number of components defined
    int rank{-1}; // rank (index) if part of a vector
    int CompFn1{0};
    int CompFn2{0};
    int isreal{1}; // T=double
    int iscomplex{0}; // T=DoubleComplex
    double CompFd1{0.0};
    double CompFd2{0.0};
    double CompFd3{0.0};
    // additional data that describe each component (defined by user):
    // occupancy, quantum number, conjugation, norm, etc.
    //Note: defined with fixed size to ease copying and MPI send
    int Nchunks[4]{0,0,0,0}; // number of chunks of each component tree
    int n1[4]{0,0,0,0};
    int n2[4]{0,0,0,0};
    int n3[4]{0,0,0,0};
    int n4[4]{0,0,0,0};
    double d1[4]{0.0,0.0,0.0,0.0};
    double d2[4]{0.0,0.0,0.0,0.0};
    double d3[4]{0.0,0.0,0.0,0.0};
};

template <int D, typename T> class CompFunction {
public:
    CompFunction(MultiResolutionAnalysis<D> &mra);
    FunctionTree<D, T> *Comp[4];

    std::string name;

    // additional data that describe each component (defined by user):
    CompFunctionData<D> data;
    int& Ncomp = data.Ncomp; //number of components defined
    int& rank = data.rank; // rank (index) if part of a vector
    int& isreal = data.isreal; // T=double
    int& iscomplex = data.iscomplex; // T=DoubleComplex
    int* Nchunks = data.Nchunks; // number of chunks of each component tree

    // ComplexFunctions are only defined for D=3
    template <int D_ = D, typename std::enable_if<D_ == 3, int>::type = 0>
    CompFunction(ComplexFunction cplxfunc);
    // CompFunction destructor
    ~CompFunction() {
        for (int i = 0; i < Ncomp; i++) {
            delete Comp[i];
        }
    }

    double norm();
    double squaredNorm();
    void alloc(int i);
    void add(T c, CompFunction<D, T> inp);
    int crop(double prec);
};

template <int D, typename T = double> using CompFunctionVector = std::vector<CompFunction<D, T> *>;

} // namespace mrcpp
