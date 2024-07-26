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
    int n1[4]{0,0,0,0}; // 0: neutral. values 1 and 2 are orthogonal to each other (product = 0)
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


template <int D> class CompFunction {
public:
    CompFunction(MultiResolutionAnalysis<D> &mra);
    CompFunction();
    CompFunction(int n1);
    CompFunction(int n1, bool share);
    CompFunction(const CompFunction<D> &compfunc);
    CompFunction(CompFunction<D> && compfunc);
    //    ComplexFunction *CPXfct; // temporary solution


    FunctionTree<D, double> *CompD[4];
    FunctionTree<D, ComplexDouble> *CompC[4];

    std::string name;

    // additional data that describe each component (defined by user):
    CompFunctionData<D> data;
    int& Ncomp = data.Ncomp; //number of components defined
    int& rank = data.rank; // rank (index) if part of a vector
    int& isreal = data.isreal; // T=double
    int& iscomplex = data.iscomplex; // T=DoubleComplex
    int* Nchunks = data.Nchunks; // number of chunks of each component tree
    // ComplexFunctions are only defined for D=3
    // template <int D_ = D, typename std::enable_if<D_ == 3, int>::type = 0>
     //CompFunction(ComplexFunction cplxfunc);
    // template <int D_ = 3, typename std::enable_if<D_ == 3, int>::type = 0>
     //operator ComplexFunction() const;
    CompFunction<D> &operator=(const CompFunction<D> &compfunc);
    // CompFunction destructor
    ~CompFunction() {
        for (int i = 0; i < Ncomp; i++) {
            delete CompD[i];
            delete CompC[i];
        }
    }

    double norm() const;
    double squaredNorm() const;
    void alloc(int i);
    void setReal(FunctionTree<D, double> *tree, int i = 0);
    void setRank(int i) {rank = i;};
    int getRank() {return rank;};
    void add(ComplexDouble c, CompFunction<D> inp);

    int crop(double prec);
    void rescale(ComplexDouble c);
    void free();
    int getSizeNodes() const;
    int getNNodes() const;

    //NB: All tbelow should be revised. Now only for backwards compatibility to ComplexFunction class
    bool hasReal()  const {return isreal;}
    bool hasImag()  const {return iscomplex;}
    bool isShared() const {return data.shared;}
    bool conjugate() const {return data.conj;}

    FunctionTree<D, double> &real() {return *CompD[0];}
    FunctionTree<D, double> &imag() {return *CompD[0];} //does not make sense
    const FunctionTree<D, double> &real() const {return *CompD[0];}
    const FunctionTree<D, double> &imag() const {return *CompD[0];} //does not make sense
    void free(int type) {delete CompD[0]; CompD[0] = nullptr; delete CompC[0]; CompC[0] = nullptr;}
    void flushFuncData();
};

template <int D> using CompFunctionVector = std::vector<CompFunction<D> *>;

template <int D>
void deep_copy(CompFunction<D> *out, const CompFunction<D> &inp);
template <int D>
void deep_copy(CompFunction<D> &out, const CompFunction<D> &inp);
template <int D>
void add(CompFunction<D> &out, ComplexDouble a, CompFunction<D> inp_a, ComplexDouble b, CompFunction<D> inp_b, double prec);
template <int D>
void linear_combination(CompFunction<D> &out, const std::vector<ComplexDouble> &c, std::vector<CompFunction<D>> &inp, double prec);
template <int D>
void multiply(CompFunction<D> &out, CompFunction<D> inp_a, CompFunction<D> inp_b, double prec, bool absPrec, bool useMaxNorms);
template <int D, typename T>
void multiply(CompFunction<D> &out, CompFunction<D> &inp_a, RepresentableFunction<D, T> &f, double prec, int nrefine = 0);
template <int D, typename T>
void multiply(CompFunction<D> &out, FunctionTree<D, T> &inp_a, RepresentableFunction<D, T> &f, double prec, int nrefine = 0);
template <int D>
ComplexDouble dot(CompFunction<D> bra, CompFunction<D> ket);
template <int D, typename T>
void project(CompFunction<D> &out, std::function<T(const Coord<D> &r)> f, double prec);
template <int D>
void project(CompFunction<D> &out, RepresentableFunction<D, double> &f, double prec);
template <int D>
void project(CompFunction<D> &out, RepresentableFunction<D, ComplexDouble> &f, double prec);

class MPI_CompFuncVector : public std::vector<CompFunction<3>> {
public:
    MPI_CompFuncVector(int N = 0);
    MultiResolutionAnalysis<3> *vecMRA;
    void distribute();
};

void rotate(MPI_CompFuncVector &Phi, const ComplexMatrix &U, double prec = -1.0);
void rotate(MPI_CompFuncVector &Phi, const ComplexMatrix &U, MPI_CompFuncVector &Psi, double prec = -1.0);
void save_nodes(MPI_CompFuncVector &Phi, mrcpp::FunctionTree<3, double> &refTree, BankAccount &account, int sizes = -1);
MPI_CompFuncVector multiply(MPI_CompFuncVector &Phi, RepresentableFunction<3> &f, double prec = -1.0, ComplexFunction *Func = nullptr, int nrefine = 1, bool all = false);
void SetdefaultMRA(MultiResolutionAnalysis<3> *MRA);
ComplexVector dot(MPI_CompFuncVector &Bra, MPI_CompFuncVector &Ket);
ComplexMatrix calc_lowdin_matrix(MPI_CompFuncVector &Phi);
ComplexMatrix calc_overlap_matrix(MPI_CompFuncVector &BraKet);
ComplexMatrix calc_overlap_matrix(MPI_CompFuncVector &Bra, MPI_CompFuncVector &Ket);
DoubleMatrix calc_norm_overlap_matrix(MPI_CompFuncVector &BraKet);
void orthogonalize(double prec, MPI_CompFuncVector &Bra, MPI_CompFuncVector &Ket);


} // namespace mrcpp
