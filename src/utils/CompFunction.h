#pragma once

#include "trees/FunctionTree.h"
#include "ComplexFunction.h"

using namespace Eigen;

namespace mrcpp {

template <int D>
struct CompFunctionData {
    // additional data that describe the overall multicomponent function (defined by user):
    // occupancy, quantum number, norm, etc.
    int Ncomp{1}; // number of components defined
    int rank{-1}; // rank (index) if part of a vector
    int conj{0}; // conjugate of all components
    int CompFn1{0};
    int CompFn2{0};
    int isreal{1}; // T=double
    int iscomplex{0}; // T=DoubleComplex
    double CompFd1{0.0};
    double CompFd2{0.0};
    double CompFd3{0.0};
    // additional data that describe each component (defined by user):
    // occupancy, quantum number, norm, etc.
    //Note: defined with fixed size to ease copying and MPI send
    int n1[4]{0,0,0,0};
    int n2[4]{0,0,0,0};
    int n3[4]{0,0,0,0};
    int n4[4]{0,0,0,0};
    double d1[4]{0.0,0.0,0.0,0.0};
    double d2[4]{0.0,0.0,0.0,0.0};
    double d3[4]{0.0,0.0,0.0,0.0};
    // used internally
    int shared{0};
    int Nchunks[4]{0,0,0,0}; // number of chunks of each component tree
};


template <int D, typename T> class CompFunction {
public:
    CompFunction();
    CompFunction(MultiResolutionAnalysis<D> &mra);
    CompFunction(const CompFunction<D, T> &compfunc);
    CompFunction(CompFunction<D, T> && compfunc);

    ComplexFunction *CPXfct; // temporary solution


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
    template <int D_ = 3, typename std::enable_if<D_ == 3, int>::type = 0>
    operator ComplexFunction() const;
    CompFunction<D, T> &operator=(const CompFunction<D, T> &func);
    // CompFunction destructor
    ~CompFunction() {
        for (int i = 0; i < Ncomp; i++) {
            delete Comp[i];
        }
    }

    double norm() const;
    double squaredNorm() const;
    void alloc(int i);
    void add(T c, CompFunction<D, T> inp);

    int crop(double prec);
    void rescale(T c);

    //NB: All tbelow should be revised. Now only for backwards compatibility to ComplexFunction class
    bool hasReal()  const {return isreal;}
    bool hasImag()  const {return iscomplex;}
    bool isShared() const {return data.shared;}

    FunctionTree<D, T> &real() {return *Comp[0];}
    FunctionTree<D, T> &imag() {return *Comp[0];}
    void free(int type) {delete Comp[0]; Comp[0] = nullptr;}
    void flushFuncData();
};

template <int D, typename T = double> using CompFunctionVector = std::vector<CompFunction<D, T> *>;

namespace compfunc {
template <int D, typename T>
void deep_copy(CompFunction<D, T> *out, const CompFunction<D, T> &inp);
template <int D, typename T>
void add(CompFunction<D, T> &out, T a, CompFunction<D, T> inp_a, T b, CompFunction<D, T> inp_b, double prec);
template <int D, typename T>
void linear_combination(CompFunction<D, T> &out, const std::vector<T> &c, std::vector<CompFunction<D, T>> &inp, double prec);
template <int D, typename T>
void multiply(CompFunction<D, T> &out, CompFunction<D, T> inp_a, CompFunction<D, T> inp_b, double prec, bool absPrec, bool useMaxNorms);
template <int D, typename T>
void multiply(CompFunction<D, T> &out, CompFunction<D, T> &inp_a, RepresentableFunction<D, T> &f, double prec, int nrefine = 0);
template <int D, typename T>
void multiply(CompFunction<D, T> &out, FunctionTree<D, T> &inp_a, RepresentableFunction<D, T> &f, double prec, int nrefine = 0);
template <int D, typename T>
T dot(CompFunction<D, T> bra, CompFunction<D, T> ket);

} // namespace compfunc
} // namespace mrcpp
