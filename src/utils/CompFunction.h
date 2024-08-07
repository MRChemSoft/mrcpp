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
    int conj{0}; // soft conjugate (all components)
    int CompFn1{0};
    int CompFn2{0};
    int isreal{0}; // trees are defined for T=double
    int iscomplex{0}; // trees are defined for T=DoubleComplex
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
    // used for storage on disk
    int type{0};
    int order{1};
    int scale{0};
    int depth{0};
    int boxes[3] = {0, 0, 0};
    int corner[3] = {0, 0, 0};

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
    CompFunction(const CompFunctionData<D>& indata);
    CompFunction(const CompFunction<D> &compfunc);
    CompFunction(CompFunction<D> && compfunc);
    CompFunction<D> &operator=(const CompFunction<D> &compfunc);

    FunctionTree<D, double> *CompD[4];
    FunctionTree<D, ComplexDouble> *CompC[4];

    std::string name;

    // additional data that describe each component (defined by user):
    CompFunctionData<D> data;
    int& Ncomp = data.Ncomp; //number of components defined
    int& rank = data.rank; // rank (index) if part of a vector
    int& conj = data.conj; // soft conjugate
    int& isreal = data.isreal; // T=double
    int& iscomplex = data.iscomplex; // T=DoubleComplex
    int* Nchunks = data.Nchunks; // number of chunks of each component tree
    // ComplexFunctions are only defined for D=3
    // template <int D_ = D, typename std::enable_if<D_ == 3, int>::type = 0>
     //CompFunction(ComplexFunction cplxfunc);
    // template <int D_ = 3, typename std::enable_if<D_ == 3, int>::type = 0>
     //operator ComplexFunction() const;
    // CompFunction destructor
    ~CompFunction() {
        for (int i = 0; i < Ncomp; i++) {
            delete CompD[i];
            delete CompC[i];
        }
    }

    CompFunction paramCopy() const;
    ComplexDouble integrate() const;
    double norm() const;
    double squaredNorm() const;
    void alloc(int i = 0);
    void setReal(FunctionTree<D, double> *tree, int i = 0);
    void setRank(int i) {rank = i;};
    const int getRank() const {return rank;};
    void add(ComplexDouble c, CompFunction<D> inp);

    int crop(double prec);
    void rescale(ComplexDouble c);
    void free();
    int getSizeNodes() const;
    int getNNodes() const;
    void flushMRAData();
    void flushFuncData();
    CompFunctionData<D> getFuncData() const;
    FunctionTree<D, double> &real(int i = 0);
    FunctionTree<D, ComplexDouble> &complex(int i = 0);
    const FunctionTree<D, double> &real(int i = 0) const;
    const FunctionTree<D, ComplexDouble> &complex(int i = 0) const;

    //NB: All below should be revised. Now only for backwards compatibility to ComplexFunction class

    void free(int type) {free();}
    bool hasReal()  const {return isreal;}
    bool hasImag()  const {return iscomplex;}
    bool isShared() const {return data.shared;}
    bool conjugate() const {return data.conj;}
    CompFunction<D> dagger();
    FunctionTree<D, double> &imag(int i = 0); //does not make sense now
    const FunctionTree<D, double> &imag(int i = 0) const; //does not make sense now
};

template <int D>
void deep_copy(CompFunction<D> *out, const CompFunction<D> &inp);
template <int D>
void deep_copy(CompFunction<D> &out, const CompFunction<D> &inp);
template <int D>
void add(CompFunction<D> &out, ComplexDouble a, CompFunction<D> inp_a, ComplexDouble b, CompFunction<D> inp_b, double prec);
template <int D>
void linear_combination(CompFunction<D> &out, const std::vector<ComplexDouble> &c, std::vector<CompFunction<D>> &inp, double prec);
template <int D>
void multiply(CompFunction<D> &out, CompFunction<D> inp_a, CompFunction<D> inp_b, double prec, bool absPrec = false, bool useMaxNorms = false);
template <int D>
void multiply(double prec, CompFunction<D> &out, double coef, CompFunction<D> inp_a, CompFunction<D> inp_b, int maxIter = -1, bool absPrec = false, bool useMaxNorms = false);
template <int D>
void multiply(CompFunction<D> &out, CompFunction<D> inp_a, CompFunction<D> inp_b, bool absPrec = false, bool useMaxNorms = false);
template <int D>
void multiply(CompFunction<D> &out, CompFunction<D> &inp_a, RepresentableFunction<D, double> &f, double prec, int nrefine = 0);
template <int D>
void multiply(CompFunction<D> &out, CompFunction<D> &inp_a, RepresentableFunction<D, ComplexDouble> &f, double prec, int nrefine = 0);
template <int D>
void multiply(CompFunction<D> &out, FunctionTree<D, double> &inp_a, RepresentableFunction<D, double> &f, double prec, int nrefine = 0);
template <int D>
void multiply(CompFunction<D> &out, FunctionTree<D, ComplexDouble> &inp_a, RepresentableFunction<D, ComplexDouble> &f, double prec, int nrefine = 0);
template <int D>
ComplexDouble dot(CompFunction<D> bra, CompFunction<D> ket);
template <int D>
double node_norm_dot(CompFunction<D> bra, CompFunction<D> ket);
void project(CompFunction<3> &out, std::function<double(const Coord<3> &r)> f, double prec);
void project(CompFunction<3> &out, std::function<ComplexDouble(const Coord<3> &r)> f, double prec);
template <int D>
void project(CompFunction<D> &out, RepresentableFunction<D, double> &f, double prec);
template <int D>
void project(CompFunction<D> &out, RepresentableFunction<D, ComplexDouble> &f, double prec);

class CompFunctionVector : public std::vector<CompFunction<3>> {
public:
    CompFunctionVector(int N = 0);
    MultiResolutionAnalysis<3> *vecMRA;
    void distribute();
};

void rotate(CompFunctionVector &Phi, const ComplexMatrix &U, double prec = -1.0);
void rotate(CompFunctionVector &Phi, const ComplexMatrix &U, CompFunctionVector &Psi, double prec = -1.0);
void save_nodes(CompFunctionVector &Phi, mrcpp::FunctionTree<3, double> &refTree, BankAccount &account, int sizes = -1);
CompFunctionVector multiply(CompFunctionVector &Phi, RepresentableFunction<3> &f, double prec = -1.0, ComplexFunction *Func = nullptr, int nrefine = 1, bool all = false);
void SetdefaultMRA(MultiResolutionAnalysis<3> *MRA);
ComplexVector dot(CompFunctionVector &Bra, CompFunctionVector &Ket);
ComplexMatrix calc_lowdin_matrix(CompFunctionVector &Phi);
ComplexMatrix calc_overlap_matrix(CompFunctionVector &BraKet);
ComplexMatrix calc_overlap_matrix(CompFunctionVector &Bra, CompFunctionVector &Ket);
DoubleMatrix calc_norm_overlap_matrix(CompFunctionVector &BraKet);
void orthogonalize(double prec, CompFunctionVector &Bra, CompFunctionVector &Ket);


} // namespace mrcpp
