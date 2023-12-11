#pragma once

#include "functions/RepresentableFunction.h"
#include "math_utils.h"
#include "mpi_utils.h"
#include "trees/FunctionTree.h"
#include "trees/MultiResolutionAnalysis.h"
#include <Eigen/Core>

using namespace Eigen;

using IntVector = Eigen::VectorXi;
using DoubleVector = Eigen::VectorXd;
using ComplexVector = Eigen::VectorXcd;

using IntMatrix = Eigen::MatrixXi;
using DoubleMatrix = Eigen::MatrixXd;
using ComplexMatrix = Eigen::MatrixXcd;

class MPI_FuncVector;

namespace mrcpp {

class BankAccount;
template <int D> class FunctionTree;
template <int D> class MultiResolutionAnalysis;

using ComplexDouble = std::complex<double>;
namespace NUMBER {
enum type { Total, Real, Imag };
}
namespace SPIN {
enum type { Paired, Alpha, Beta };
}

struct FunctionData {
    int type{0};
    int order{1};
    int scale{0};
    int depth{0};
    int boxes[3] = {0, 0, 0};
    int corner[3] = {0, 0, 0};
    int real_size{0};
    int imag_size{0};
    bool is_shared{false};
    int spin{0};
    int occ{0};
};

class TreePtr final {
public:
    explicit TreePtr(bool share)
            : shared_mem_re(nullptr)
            , shared_mem_im(nullptr)
            , re(nullptr)
            , im(nullptr) {
        this->func_data.is_shared = share;
        if (this->func_data.is_shared and mpi::share_size > 1) {
            // Memory size in MB defined in input. Virtual memory, does not cost anything if not used.
#ifdef MRCPP_HAS_MPI
            this->shared_mem_re = new mrcpp::SharedMemory(mpi::comm_share, mpi::shared_memory_size);
            this->shared_mem_im = new mrcpp::SharedMemory(mpi::comm_share, mpi::shared_memory_size);
#endif
        }
    }

    ~TreePtr() {
        if (this->shared_mem_re != nullptr) delete this->shared_mem_re;
        if (this->shared_mem_im != nullptr) delete this->shared_mem_im;
        if (this->re != nullptr) delete this->re;
        if (this->im != nullptr) delete this->im;
    }

    friend class ComplexFunction;

private:
    FunctionData func_data;
    mrcpp::SharedMemory *shared_mem_re;
    mrcpp::SharedMemory *shared_mem_im;
    mrcpp::FunctionTree<3> *re; ///< Real part of function
    mrcpp::FunctionTree<3> *im; ///< Imaginary part of function

    void flushFuncData() {
        this->func_data.real_size = 0;
        this->func_data.imag_size = 0;
        if (this->re != nullptr) {
            this->func_data.real_size = this->re->getNChunksUsed();
            flushMRAData(this->re->getMRA());
        }
        if (this->im != nullptr) {
            this->func_data.imag_size = this->im->getNChunksUsed();
            flushMRAData(this->im->getMRA());
        }
    }

    void flushMRAData(const mrcpp::MultiResolutionAnalysis<3> &mra) {
        const auto &box = mra.getWorldBox();
        this->func_data.type = mra.getScalingBasis().getScalingType();
        this->func_data.order = mra.getOrder();
        this->func_data.depth = mra.getMaxDepth();
        this->func_data.scale = box.getScale();
        this->func_data.boxes[0] = box.size(0);
        this->func_data.boxes[1] = box.size(1);
        this->func_data.boxes[2] = box.size(2);
        this->func_data.corner[0] = box.getCornerIndex().getTranslation(0);
        this->func_data.corner[1] = box.getCornerIndex().getTranslation(1);
        this->func_data.corner[2] = box.getCornerIndex().getTranslation(2);
    }
};

class ComplexFunction {
public:
    ComplexFunction(std::shared_ptr<TreePtr> funcptr);
    ComplexFunction(const ComplexFunction &func);
    ComplexFunction(int spin = 0, int occ = -1, int rank = -1, bool share = false);
    ComplexFunction &operator=(const ComplexFunction &func);
    ComplexFunction paramCopy() const;
    bool isShared() const { return this->func_ptr->func_data.is_shared; }
    bool hasReal() const { return (this->func_ptr->re == nullptr) ? false : true; }
    bool hasImag() const { return (this->func_ptr->im == nullptr) ? false : true; }
    FunctionData &getFunctionData();
    int occ() const { return this->func_ptr->func_data.occ; }
    int spin() const { return this->func_ptr->func_data.spin; }
    FunctionTree<3> &real() { return *this->func_ptr->re; }
    FunctionTree<3> &imag() { return *this->func_ptr->im; }
    const FunctionTree<3> &real() const { return *this->func_ptr->re; }
    const FunctionTree<3> &imag() const { return *this->func_ptr->im; }
    void release() { this->func_ptr.reset(); }
    bool conjugate() const { return this->conj; }
    MultiResolutionAnalysis<3> *funcMRA = nullptr;
    int getRank() const { return rank; }
    void setRank(int rank) { (*this).rank = rank; }
    void setOcc(int occ) { this->getFunctionData().occ = occ; }
    void setSpin(int spin) { this->getFunctionData().spin = spin; }
    ComplexFunction dagger();
    virtual ~ComplexFunction() = default;

    void alloc(int type, mrcpp::MultiResolutionAnalysis<3> *mra = nullptr);
    void free(int type);

    int getSizeNodes(int type) const;
    int getNNodes(int type) const;

    void setReal(mrcpp::FunctionTree<3> *tree);
    void setImag(mrcpp::FunctionTree<3> *tree);

    double norm() const;
    double squaredNorm() const;
    ComplexDouble integrate() const;

    int crop(double prec);
    void rescale(double c);
    void rescale(ComplexDouble c);
    void add(ComplexDouble c, ComplexFunction inp);
    void absadd(ComplexDouble c, ComplexFunction inp);
    char printSpin() const;

protected:
    bool conj{false};
    std::shared_ptr<mrcpp::TreePtr> func_ptr;
    int rank = -1; // index in vector
};

namespace cplxfunc {
void SetdefaultMRA(MultiResolutionAnalysis<3> *MRA);
ComplexDouble dot(ComplexFunction bra, ComplexFunction ket);
ComplexDouble node_norm_dot(ComplexFunction bra, ComplexFunction ket, bool exact);
void deep_copy(ComplexFunction &out, ComplexFunction &inp);
void add(ComplexFunction &out, ComplexDouble a, ComplexFunction inp_a, ComplexDouble b, ComplexFunction inp_b, double prec);
void project(ComplexFunction &out, std::function<double(const Coord<3> &r)> f, int type, double prec);
void project(ComplexFunction &out, RepresentableFunction<3> &f, int type, double prec);
void multiply(ComplexFunction &out, ComplexFunction inp_a, ComplexFunction inp_b, double prec, bool absPrec = false, bool useMaxNorms = false);
void multiply_real(ComplexFunction &out, ComplexFunction inp_a, ComplexFunction inp_b, double prec, bool absPrec = false, bool useMaxNorms = false);
void multiply_imag(ComplexFunction &out, ComplexFunction inp_a, ComplexFunction inp_b, double prec, bool absPrec = false, bool useMaxNorms = false);
void multiply(ComplexFunction &out, ComplexFunction &inp_a, RepresentableFunction<3> &f, double prec, int nrefine = 0);
void multiply(ComplexFunction &out, FunctionTree<3> &inp_a, RepresentableFunction<3> &f, double prec, int nrefine = 0);
void linear_combination(ComplexFunction &out, const ComplexVector &c, std::vector<ComplexFunction> &inp, double prec);
} // namespace cplxfunc

class MPI_FuncVector : public std::vector<ComplexFunction> {
public:
    MPI_FuncVector(int N = 0);
    MultiResolutionAnalysis<3> *vecMRA;
    void distribute();
};

namespace mpifuncvec {
void rotate(MPI_FuncVector &Phi, const ComplexMatrix &U, double prec = -1.0);
void rotate(MPI_FuncVector &Phi, const ComplexMatrix &U, MPI_FuncVector &Psi, double prec = -1.0);
void save_nodes(MPI_FuncVector &Phi, mrcpp::FunctionTree<3> &refTree, BankAccount &account, int sizes = -1);
MPI_FuncVector multiply(MPI_FuncVector &Phi, RepresentableFunction<3> &f, double prec = -1.0, ComplexFunction *Func = nullptr, int nrefine = 1, bool all = false);
ComplexVector dot(MPI_FuncVector &Bra, MPI_FuncVector &Ket);
ComplexMatrix calc_lowdin_matrix(MPI_FuncVector &Phi);
ComplexMatrix calc_overlap_matrix(MPI_FuncVector &BraKet);
ComplexMatrix calc_overlap_matrix(MPI_FuncVector &Bra, MPI_FuncVector &Ket);
DoubleMatrix calc_norm_overlap_matrix(MPI_FuncVector &BraKet);
void orthogonalize(double prec, MPI_FuncVector &Bra, MPI_FuncVector &Ket);
} // namespace mpifuncvec
} // namespace mrcpp
