#include "CompFunction.h"
#include "Bank.h"
#include "Printer.h"
#include "parallel.h"
#include "treebuilders/add.h"
#include "treebuilders/grid.h"
#include "treebuilders/multiply.h"
#include "treebuilders/project.h"
#include "trees/FunctionNode.h"
#include <fstream>

/* Some rules for CompFunction:
 * NComp is the number of components. If Ncomp>0, the corresponding trees must exist (can be only empty roots).
 * The other trees should be set to nullptr.
 * The trees and data can be shared among several CompFunction; this is managed automatically by "std::make_shared"
 * Normally the CompFunction must be eiher real or complex (or none if noe is defined anyway).
 * Though it is allowed in some cases to have both and the code should preferably allow this. (It is used temporary
 * when we need a Complex type, but the trees are real: the tree is then copied as a complex tree in the same CompFunction).
 * TreePtr (aka func_ptr) is the part potentially shared with others with "std::make_shared". It contains the pointers to the trees.
 * The static data (number of components, real/complex, conjugaison, integers used for spin etc.) are store in func_ptr.data.
 */

namespace mrcpp {

template <int D> MultiResolutionAnalysis<D> *defaultCompMRA = nullptr; // Global MRA

template <int D> CompFunction<D>::CompFunction(MultiResolutionAnalysis<D> &mra) {
    defaultCompMRA<D> = &mra;
    func_ptr = std::make_shared<TreePtr<D>>(false);
    CompD = func_ptr->real;
    CompC = func_ptr->cplx;
    for (int i = 0; i < 4; i++) CompD[i] = nullptr;
    for (int i = 0; i < 4; i++) CompC[i] = nullptr;
}

template <int D> CompFunction<D>::CompFunction() {
    func_ptr = std::make_shared<TreePtr<D>>(false);
    CompD = func_ptr->real;
    CompC = func_ptr->cplx;
    for (int i = 0; i < 4; i++) CompD[i] = nullptr;
    for (int i = 0; i < 4; i++) CompC[i] = nullptr;
}

/*
 * Empty functions (no components defined)
 */
template <int D> CompFunction<D>::CompFunction(int n1) {
    func_ptr = std::make_shared<TreePtr<D>>(false);
    CompD = func_ptr->real;
    CompC = func_ptr->cplx;
    for (int i = 0; i < 4; i++) CompD[i] = nullptr;
    for (int i = 0; i < 4; i++) CompC[i] = nullptr;
    func_ptr->data.n1[0] = n1;
    func_ptr->data.n2[0] = -1;
    func_ptr->data.n3[0] = 0;
    func_ptr->rank = 0;
    func_ptr->isreal = 1;
    func_ptr->iscomplex = 0;
    func_ptr->data.shared = false;
}

/*
 * Empty functions (no components defined)
 */
template <int D> CompFunction<D>::CompFunction(int n1, bool share) {
    func_ptr = std::make_shared<TreePtr<D>>(share);
    CompD = func_ptr->real;
    CompC = func_ptr->cplx;
    for (int i = 0; i < 4; i++) CompD[i] = nullptr;
    for (int i = 0; i < 4; i++) CompC[i] = nullptr;
    func_ptr->data.n1[0] = n1;
    func_ptr->data.n2[0] = -1;
    func_ptr->data.n3[0] = 0;
    func_ptr->rank = 0;
    func_ptr->isreal = 1;
    func_ptr->iscomplex = 0;
    func_ptr->data.shared = share;
}

/*
 * Empty functions (trees defined but zero)
 */
template <int D> CompFunction<D>::CompFunction(const CompFunctionData<D> &indata, bool alloc) {
    func_ptr = std::make_shared<TreePtr<D>>(indata.shared);
    func_ptr->data = indata;
    CompD = func_ptr->real;
    CompC = func_ptr->cplx;
    if (alloc)
        this->alloc(Ncomp());
    else
        this->free();
}

/** @brief Copy constructor
 *
 * Shallow copy: meta data is copied along with the component pointers,
 * NO transfer of ownership.
 */
template <int D> CompFunction<D>::CompFunction(const CompFunction<D> &compfunc) {
    func_ptr = compfunc.func_ptr;
    CompD = func_ptr->real;
    CompC = func_ptr->cplx;
}

/** @brief Copy constructor
 *
 * Shallow copy: meta data is copied along with the component pointers,
 * NO transfer of ownership.
 */
template <int D> CompFunction<D> &CompFunction<D>::operator=(const CompFunction<D> &compfunc) {
    if (this != &compfunc) {
        func_ptr = compfunc.func_ptr;
        CompD = func_ptr->real;
        CompC = func_ptr->cplx;
    }
    return *this;
}

template <int D>
/** @brief Parameter copy
 *
 * Returns a copy without defined trees.
 */
CompFunction<D> CompFunction<D>::paramCopy(bool alloc) const {
    CompFunction<D> out(func_ptr->data, alloc);
    // we do not copy tree sizes:
    for (int i = 0; i < 4; i++) out.func_ptr->data.Nchunks[i] = 0;
    return out;
}

template <int D> void CompFunction<D>::flushMRAData() {
    const auto &box = defaultCompMRA<3>->getWorldBox();
    func_ptr->data.type = defaultCompMRA<3>->getScalingBasis().getScalingType();
    func_ptr->data.order = defaultCompMRA<3>->getOrder();
    func_ptr->data.depth = defaultCompMRA<3>->getMaxDepth();
    func_ptr->data.scale = box.getScale();
    func_ptr->data.boxes[0] = box.size(0);
    func_ptr->data.boxes[1] = box.size(1);
    func_ptr->data.boxes[2] = box.size(2);
    func_ptr->data.corner[0] = box.getCornerIndex().getTranslation(0);
    func_ptr->data.corner[1] = box.getCornerIndex().getTranslation(1);
    func_ptr->data.corner[2] = box.getCornerIndex().getTranslation(2);
}

template <int D> void CompFunction<D>::flushFuncData() {
    if (D == 3) flushMRAData();
    for (int i = 0; i < Ncomp(); i++) {
        if (isreal()) {
            func_ptr->Nchunks[i] = CompD[i]->getNChunksUsed();
        } else {
            func_ptr->Nchunks[i] = CompC[i]->getNChunksUsed();
        }
    }
    for (int i = Ncomp(); i < 4; i++) func_ptr->Nchunks[i] = 0;
}

template <int D> CompFunctionData<D> CompFunction<D>::getFuncData() const {
    CompFunctionData<D> outdata;
    const auto &box = defaultCompMRA<3>->getWorldBox();
    outdata.type = defaultCompMRA<3>->getScalingBasis().getScalingType();
    outdata.order = defaultCompMRA<3>->getOrder();
    outdata.depth = defaultCompMRA<3>->getMaxDepth();
    outdata.scale = box.getScale();
    outdata.boxes[0] = box.size(0);
    outdata.boxes[1] = box.size(1);
    outdata.boxes[2] = box.size(2);
    outdata.corner[0] = box.getCornerIndex().getTranslation(0);
    outdata.corner[1] = box.getCornerIndex().getTranslation(1);
    outdata.corner[2] = box.getCornerIndex().getTranslation(2);
    for (int i = 0; i < Ncomp(); i++) {
        if (isreal()) {
            outdata.Nchunks[i] = CompD[i]->getNChunksUsed();
        } else {
            outdata.Nchunks[i] = CompC[i]->getNChunksUsed();
        }
    }
    for (int i = Ncomp(); i < 4; i++) outdata.Nchunks[i] = 0;
    return outdata;
}

template <int D> ComplexDouble CompFunction<D>::integrate() const {
    ComplexDouble integral;
    if (isreal())
        integral = CompD[0]->integrate();
    else
        integral = CompC[0]->integrate();
    return integral;
}

template <int D> double CompFunction<D>::norm() const {
    double norm = getSquareNorm();
    if (norm > 0.0) norm = std::sqrt(norm);
    return norm;
}
template <int D> double CompFunction<D>::getSquareNorm() const {
    double norm = 0.0;
    for (int i = 0; i < Ncomp(); i++) {
        if (isreal() and CompD[i] != nullptr) {
            norm += CompD[i]->getSquareNorm();
        } else if (iscomplex() and CompC[i] != nullptr) {
            norm += CompC[i]->getSquareNorm();
        }
    }
    return norm;
}

//  Allocate empty trees. The tree must be defined as real or complex already.
//  Allocates all ialloc trees, with indices 0,...ialloc-1
//  nalloc is the number of components allocated. ialloc=1 allocates one tree.
//  deletes all old trees if found.
template <int D> void CompFunction<D>::alloc(int nalloc, bool zero) {
    if (defaultCompMRA<D> == nullptr) MSG_ABORT("Default MRA not yet defined");
    if (isreal() == 0 and iscomplex() == 0) MSG_ABORT("Function must be defined either real or complex");
    for (int i = 0; i < nalloc; i++) {
        delete CompD[i];
        delete CompC[i];
        CompD[i] = nullptr;
        CompC[i] = nullptr;
        if (isreal()) {
            CompD[i] = new FunctionTree<D, double>(*defaultCompMRA<D>, func_ptr->shared_mem_real);
            if (zero) CompD[i]->setZero();
        }
        if (iscomplex()) {
            CompC[i] = new FunctionTree<D, ComplexDouble>(*defaultCompMRA<D>, func_ptr->shared_mem_cplx);
            if (zero) CompC[i]->setZero();
        }
        func_ptr->Ncomp = std::max(Ncomp(), i + 1);
    }
    for (int i = nalloc; i < Ncomp(); i++) {
        // delete possible remaining components
        delete CompD[i];
        delete CompC[i];
        CompD[i] = nullptr;
        CompC[i] = nullptr;
    }
}

//  Allocate one empty trees for one specific component.
//  The tree must be defined as real or complex already.
//  ialloc is index allocated. ialloc=0 allocates the tree with index zero.
//  deletes old tree if found.
template <int D> void CompFunction<D>::alloc_comp(int ialloc) {
    if (defaultCompMRA<D> == nullptr) MSG_ABORT("Default MRA not yet defined");
    if (isreal() == 0 and iscomplex() == 0) MSG_ABORT("Function must be defined either real or complex");
    int i = ialloc;
    delete CompD[i];
    delete CompC[i];
    CompD[i] = nullptr;
    CompC[i] = nullptr;
    if (isreal()) {
        CompD[i] = new FunctionTree<D, double>(*defaultCompMRA<D>, func_ptr->shared_mem_real);
        CompD[i]->setZero();
    }
    if (iscomplex()) {
        CompC[i] = new FunctionTree<D, ComplexDouble>(*defaultCompMRA<D>, func_ptr->shared_mem_cplx);
        CompC[i]->setZero();
    }
    func_ptr->Ncomp = std::max(Ncomp(), i + 1);
}

template <int D> void CompFunction<D>::free() {
    for (int i = 0; i < Ncomp(); i++) {
        if (CompD[i] != nullptr) delete CompD[i];
        if (CompC[i] != nullptr) delete CompC[i];
        CompD[i] = nullptr;
        CompC[i] = nullptr;
    }
    if (this->func_ptr->shared_mem_real) this->func_ptr->shared_mem_real->clear();
    if (this->func_ptr->shared_mem_cplx) this->func_ptr->shared_mem_cplx->clear();
    func_ptr->Ncomp = 0;
}

template <int D> int CompFunction<D>::getSizeNodes() const {
    int size_mb = 0; // Memory size in kB
    for (int i = 0; i < Ncomp(); i++) {
        if (isreal() and CompD[i] != nullptr) size_mb += CompD[i]->getSizeNodes();
        if (iscomplex() and CompC[i] != nullptr) size_mb += CompC[i]->getSizeNodes();
    }
    return size_mb;
}

template <int D> int CompFunction<D>::getNNodes() const {
    int nNodes = 0;
    for (int i = 0; i < Ncomp(); i++) {
        if (isreal() and CompD[i] != nullptr) nNodes += CompD[i]->getNNodes();
        if (iscomplex() and CompC[i] != nullptr) nNodes += CompC[i]->getNNodes();
    }
    return nNodes;
}

/** @brief Soft complex conjugate
 *
 * Will use complex conjugate in operations (add, multiply etc.)
 * Does change the state (conj flag), but does not actively change all coefficients.
 */
template <int D> void CompFunction<D>::dagger() {
    func_ptr->data.conj = not(func_ptr->data.conj);
    for (int i = 0; i < Ncomp(); i++) {
        if (CompC[i] != nullptr) CompC[i]->setConjugate(func_ptr->data.conj);
    }
}

template <int D> FunctionTree<D, double> &CompFunction<D>::real(int i) {
    if (!isreal()) MSG_ABORT("not real function");
    if (CompD[i] == nullptr) alloc_comp(i);
    return *CompD[i];
}
template <int D> // NB: should return CompC in the future
FunctionTree<D, double> &CompFunction<D>::imag(int i) {
    MSG_ABORT("Must choose real or complex");
    if (!iscomplex()) MSG_ABORT("not complex function");
    return *CompD[i];
}

template <int D> FunctionTree<D, ComplexDouble> &CompFunction<D>::complex(int i) {
    if (!iscomplex()) MSG_ABORT("not marked as a complex function");
    if (CompC[i] == nullptr) alloc_comp(i);
    return *CompC[i];
}

template <int D> const FunctionTree<D, double> &CompFunction<D>::real(int i) const {
    if (!isreal()) MSG_ABORT("not real function");
    return *CompD[i];
}
template <int D> // NB: should use complex or real
const FunctionTree<D, double> &CompFunction<D>::imag(int i) const {
    MSG_ABORT("Must choose real or complex");
    if (!iscomplex()) MSG_ABORT("not complex function");
    return *CompD[i];
}
template <int D> const FunctionTree<D, ComplexDouble> &CompFunction<D>::complex(int i) const {
    if (!iscomplex()) MSG_ABORT("not marked as a complex function");
    return *CompC[i];
}

/* for backwards compatibility */
template <int D> void CompFunction<D>::setReal(FunctionTree<D, double> *tree, int i) {
    func_ptr->isreal = 1;
    // if (CompD[i] != nullptr) delete CompD[i];
    CompD[i] = tree;
    if (tree != nullptr) {
        func_ptr->Ncomp = std::max(Ncomp(), i + 1);
    } else {
        func_ptr->Ncomp = std::min(Ncomp(), i);
    }
}

template <int D> void CompFunction<D>::setCplx(FunctionTree<D, ComplexDouble> *tree, int i) {
    func_ptr->iscomplex = 1;
    // if (CompC[i] != nullptr) delete CompC[i];
    CompC[i] = tree;
    if (tree != nullptr) {
        func_ptr->Ncomp = std::max(Ncomp(), i + 1);
    } else {
        func_ptr->Ncomp = std::min(Ncomp(), i);
    }
}

/** @brief In place addition.
 *
 * Output is extended to union grid.
 *
 */
template <int D> void CompFunction<D>::add(ComplexDouble c, CompFunction<D> inp) {

    if (Ncomp() < inp.Ncomp()) {
        func_ptr->data = inp.func_ptr->data;
        alloc(inp.Ncomp(), true);
    }

    for (int i = 0; i < inp.Ncomp(); i++) {
        if (inp.isreal() and c.imag() < MachineZero) {
            CompD[i]->add_inplace(c.real(), *inp.CompD[i]);
        } else {
            if (this->isreal()) {
                CompD[i]->CopyTreeToComplex(CompC[i]);
                delete CompD[i];
                CompD[i] = nullptr;
                func_ptr->iscomplex = true;
                func_ptr->isreal = false;
            }
            CompC[i]->add_inplace(c, *inp.CompC[i]);
        }
    }
}

template <int D> int CompFunction<D>::crop(double prec) {
    if (prec < 0.0) return 0;
    int nChunksremoved = 0;
    for (int i = 0; i < Ncomp(); i++) {
        if (isreal()) {
            nChunksremoved += CompD[i]->crop(prec, 1.0, false);
        } else {
            nChunksremoved += CompC[i]->crop(prec, 1.0, false);
        }
    }
    return nChunksremoved;
}

/** @brief In place multiply with scalar. Fully in-place.*/
template <int D> void CompFunction<D>::rescale(ComplexDouble c) {
    bool need_to_rescale = not(isShared()) or mpi::share_master();
    if (need_to_rescale) {
        for (int i = 0; i < Ncomp(); i++) {
            if (iscomplex()) {
                CompC[i]->rescale(c);
            } else {
                if (abs(c.imag()) > MachineZero) { // works only only for NComp==1)
                    CompD[i]->CopyTreeToComplex(CompC[i]);
                    delete CompD[i];
                    CompD[i] = nullptr;
                    func_ptr->iscomplex = true;
                    func_ptr->isreal = false;
                    CompC[i]->rescale(c);
                } else {
                    CompD[i]->rescale(c.real());
                }
            }
        }
    } else
        MSG_ABORT("Not implemented");
}

template class MultiResolutionAnalysis<1>;
template class MultiResolutionAnalysis<2>;
template class MultiResolutionAnalysis<3>;
template class CompFunction<1>;
template class CompFunction<2>;
template class CompFunction<3>;

/** @brief Deep copy
 *
 * Deep copy: meta data is copied along with the content of each component.
 */
template <int D> void deep_copy(CompFunction<D> *out, const CompFunction<D> &inp) {
    out->func_ptr->data = inp.func_ptr->data;
    out->alloc(inp.Ncomp());
    if (inp.getNNodes() == 0) return;
    for (int i = 0; i < inp.Ncomp(); i++) {
        if (inp.isreal()) {
            inp.CompD[i]->deep_copy(out->CompD[i]);
        } else {
            inp.CompC[i]->deep_copy(out->CompC[i]);
        }
    }
}

/** @brief Deep copy
 *
 * Deep copy: meta func_ptr->data is copied along with the content of each component.
 */
template <int D> void deep_copy(CompFunction<D> &out, const CompFunction<D> &inp) {
    out.func_ptr->data = inp.func_ptr->data;
    out.alloc(inp.Ncomp());
    if (inp.getNNodes() == 0) return;
    for (int i = 0; i < inp.Ncomp(); i++) {
        if (inp.isreal()) {
            inp.CompD[i]->deep_copy(out.CompD[i]);
        } else {
            inp.CompC[i]->deep_copy(out.CompC[i]);
        }
    }
}

/** @brief out = a*inp_a + b*inp_b
 *
 * Recast into linear_combination.
 *
 */
template <int D> void add(CompFunction<D> &out, ComplexDouble a, CompFunction<D> inp_a, ComplexDouble b, CompFunction<D> inp_b, double prec, bool conjugate) {
    std::vector<ComplexDouble> coefs(2);
    coefs[0] = a;
    coefs[1] = b;

    std::vector<CompFunction<D>> funcs; // NB: not a CompFunctionVector, because not run in parallel!
    funcs.push_back(inp_a);
    funcs.push_back(inp_b);

    linear_combination(out, coefs, funcs, prec, conjugate);
}

/** @brief out = c_0*inp_0 + c_1*inp_1 + ... + c_N*inp_N
 *
 * OMP parallel, but not MPI parallel
 */
template <int D> void linear_combination(CompFunction<D> &out, const std::vector<ComplexDouble> &c, std::vector<CompFunction<D>> &inp, double prec, bool conjugate) {
    double thrs = MachineZero;
    bool need_to_add = not(out.isShared()) or mpi::share_master();
    bool share = out.isShared();
    out.func_ptr->data = inp[0].func_ptr->data;
    out.func_ptr->data.shared = share; // we don' inherit the shareness
    bool iscomplex = false;
    for (int i = 0; i < inp.size(); i++)
        if (inp[i].iscomplex() or c[i].imag() > MachineZero) iscomplex = true;
    if (iscomplex) {
        out.func_ptr->data.iscomplex = 1;
        out.func_ptr->data.isreal = 0;
    }
    out.alloc(out.Ncomp());
    for (int comp = 0; comp < inp[0].Ncomp(); comp++) {
        if (not iscomplex) {
            FunctionTreeVector<D, double> fvec; // one component vector
            for (int i = 0; i < inp.size(); i++) {
                if (std::norm(c[i]) < thrs) continue;
                if (inp[i].getNNodes() == 0 or inp[i].CompD[comp]->getSquareNorm() < thrs) continue;
                fvec.push_back(std::make_tuple(c[i].real(), inp[i].CompD[comp]));
            }
            if (need_to_add) {
                if (fvec.size() > 0) {
                    if (prec < 0.0) {
                        build_grid(*out.CompD[comp], fvec);
                        mrcpp::add(prec, *out.CompD[comp], fvec, 0);
                    } else {
                        mrcpp::add(prec, *out.CompD[comp], fvec);
                    }
                } else if (out.isreal()) {
                    out.CompD[comp]->setZero();
                }
            }
        } else {
            FunctionTreeVector<D, ComplexDouble> fvec; // one component vector
            for (int i = 0; i < inp.size(); i++) {
                if (inp[i].isreal()) {
                    inp[i].CompD[comp]->CopyTreeToComplex(inp[i].CompC[comp]);
                    delete inp[i].CompD[comp];
                    inp[i].CompD[comp] = nullptr;
                    inp[i].func_ptr->iscomplex = true;
                    inp[i].func_ptr->isreal = false;
                }
                if (std::norm(c[i]) < thrs) continue;
                if (inp[i].getNNodes() == 0 or inp[i].CompC[comp]->getSquareNorm() < thrs) continue;
                fvec.push_back(std::make_tuple(c[i], inp[i].CompC[comp]));
            }
            if (need_to_add) {
                if (fvec.size() > 0) {
                    if (prec < 0.0) {
                        build_grid(*out.CompC[comp], fvec);
                        mrcpp::add(prec, *out.CompC[comp], fvec, 0, false, conjugate);
                    } else {
                        mrcpp::add(prec, *out.CompC[comp], fvec, -1, false, conjugate);
                    }
                } else if (out.iscomplex()) {
                    out.CompC[comp]->setZero();
                }
            }
        }
        mpi::share_function(out, 0, 9911, mpi::comm_share);
    }
}

/** @brief out = inp_a * inp_b
 *
 */
template <int D> void multiply(CompFunction<D> &out, CompFunction<D> inp_a, CompFunction<D> inp_b, double prec, bool absPrec, bool useMaxNorms, bool conjugate) {
    multiply(prec, out, 1.0, inp_a, inp_b, -1, absPrec, useMaxNorms, conjugate);
}

/** @brief out = inp_a * inp_b
 *  Takes conjugate of inp_a if conjugate=true
 *  In case of mixed real/complex inputs, the real functions are converted into complex functions.
 */
template <int D> void multiply(double prec, CompFunction<D> &out, double coef, CompFunction<D> inp_a, CompFunction<D> inp_b, int maxIter, bool absPrec, bool useMaxNorms, bool conjugate) {
    if (inp_b.func_ptr->conj) MSG_ABORT("Not implemented");
    if (inp_a.func_ptr->conj) conjugate = (not conjugate);
    bool need_to_multiply = not(out.isShared()) or mpi::share_master();
    bool out_allocated = true;
    if (out.Ncomp() == 0) out_allocated = false;
    bool share = out.isShared();
    out.func_ptr->data = inp_a.func_ptr->data;
    out.func_ptr->data.shared = share; // we don't inherit the shareness
    out.func_ptr->conj = false;        // we don't inherit conjugaison
    if (inp_a.getNNodes() == 0 or inp_b.getNNodes() == 0) {
        if (!out_allocated) out.alloc(out.Ncomp());
        return;
    }
    for (int comp = 0; comp < inp_a.Ncomp(); comp++) {
        out.func_ptr->data.c1[comp] = inp_a.func_ptr->data.c1[comp] * inp_b.func_ptr->data.c1[comp]; // we could put this is coef if everything is real?
        if (inp_a.isreal() and inp_b.isreal()) {
            if (need_to_multiply) {
                if (!out_allocated) out.alloc(out.Ncomp());
                if (prec < 0.0) {
                    // Union grid
                    build_grid(*out.CompD[comp], *inp_a.CompD[comp]);
                    build_grid(*out.CompD[comp], *inp_b.CompD[comp]);
                    mrcpp::multiply(prec, *out.CompD[comp], coef, *inp_a.CompD[comp], *inp_b.CompD[comp], 0, false, false, conjugate);
                } else {
                    // Adaptive grid
                    mrcpp::multiply(prec, *out.CompD[comp], coef, *inp_a.CompD[comp], *inp_b.CompD[comp], maxIter, absPrec, useMaxNorms, conjugate);
                }
            }
        } else {
            // if one of the input is real, we simply make a new complex copy of it
            bool inp_aisReal = inp_a.isreal();
            bool inp_bisReal = inp_b.isreal();
            if (inp_aisReal) {
                inp_a.CompD[comp]->CopyTreeToComplex(inp_a.CompC[comp]);
                inp_a.func_ptr->iscomplex = true;
                inp_a.func_ptr->isreal = false;
            }
            if (inp_bisReal) {
                inp_b.CompD[comp]->CopyTreeToComplex(inp_b.CompC[comp]);
                inp_b.func_ptr->iscomplex = true;
                inp_b.func_ptr->isreal = false;
            }
            ComplexDouble coef = 1.0;
            if (need_to_multiply) {
                if (prec < 0.0) {
                    // Union grid
                    out.func_ptr->iscomplex = 1;
                    out.func_ptr->isreal = 0;
                    delete out.CompD[comp];
                    delete out.CompC[comp];
                    if (!out_allocated) out.alloc(out.Ncomp());
                    build_grid(*out.CompC[comp], *inp_a.CompC[comp]);
                    build_grid(*out.CompC[comp], *inp_b.CompC[comp]);
                    mrcpp::multiply(prec, *out.CompC[comp], coef, *inp_a.CompC[comp], *inp_b.CompC[comp], 0, false, false, conjugate);
                } else { // note that this assumes Ncomp=1
                    // Adaptive grid
                    if (out.CompD[comp] != nullptr) { // NB: func_ptr has alreadybeen overwritten!
                        if (out.CompD[comp]->getNNodes() > 0) {
                            out.CompD[comp]->CopyTreeToComplex(out.CompC[comp]);
                            out.func_ptr->iscomplex = 1;
                            out.func_ptr->isreal = 0;
                            delete out.CompD[comp];
                            out.CompD[comp] = nullptr;
                        } else {
                            out.func_ptr->iscomplex = 1;
                            out.func_ptr->isreal = 0;
                            out.alloc(out.Ncomp());
                        }
                    } else {
                        out.func_ptr->iscomplex = 1;
                        out.func_ptr->isreal = 0;
                        if (!out_allocated) out.alloc(out.Ncomp());
                    }
                    mrcpp::multiply(prec, *out.CompC[comp], coef, *inp_a.CompC[comp], *inp_b.CompC[comp], maxIter, absPrec, useMaxNorms, conjugate);
                }
            }
            // restore original tree
            if (inp_aisReal) {
                delete inp_a.CompC[comp];
                inp_a.CompC[comp] = nullptr;
                inp_a.func_ptr->iscomplex = false;
                inp_a.func_ptr->isreal = true;
            }
            if (inp_bisReal) {
                delete inp_b.CompC[comp];
                inp_b.CompC[comp] = nullptr;
                inp_b.func_ptr->iscomplex = false;
                inp_b.func_ptr->isreal = true;
            }
        }
    }
    mpi::share_function(out, 0, 9911, mpi::comm_share);
}

/** @brief out = inp_a * f
 *
 *  Only one component is multiplied
 */
template <int D> void multiply(CompFunction<D> &out, CompFunction<D> &inp_a, RepresentableFunction<D, double> &f, double prec, int nrefine, bool conjugate) {
    if (inp_a.Ncomp() > 1) MSG_ABORT("Not implemented");
    if (inp_a.isreal() != 1) MSG_ABORT("Not implemented");
    if (conjugate) MSG_ABORT("Not implemented");
    CompFunctionVector CompVec; // Should use vector<CompFunction>?
    CompVec.push_back(inp_a);
    CompFunctionVector CompVecOut;
    CompVecOut = multiply(CompVec, f, prec, nullptr, nrefine, true);
    out = CompVecOut[0];
    //    multiply(out, *inp_a.CompD[0], f, prec, nrefine, conjugate);
}

/** @brief out = inp_a * f
 *
 *  Only one component is multiplied
 */
template <int D> void multiply(CompFunction<D> &out, CompFunction<D> &inp_a, RepresentableFunction<D, ComplexDouble> &f, double prec, int nrefine, bool conjugate) {
    MSG_ABORT("Not implemented");
    if (inp_a.Ncomp() > 1) MSG_ABORT("Not implemented");
    if (inp_a.iscomplex() != 1) MSG_ABORT("Not implemented");
    if (conjugate) MSG_ABORT("Not implemented");
    CompFunctionVector CompVec; // Should use vector<CompFunction>?
    CompVec.push_back(inp_a);
    CompFunctionVector CompVecOut;
    // CompVecOut = multiply(CompVec, f, prec, nrefine, true);
    out = CompVecOut[0];
}

/** @brief out = inp_a * f
 *
 */
template <int D> void multiply(CompFunction<D> &out, FunctionTree<D, double> &inp_a, RepresentableFunction<D, double> &f, double prec, int nrefine, bool conjugate) {
    CompFunction<D> func_a;
    func_a.func_ptr->isreal = 1;
    func_a.func_ptr->iscomplex = 0;
    func_a.alloc(1);
    func_a.CompD[0] = &inp_a;
    multiply(out, func_a, f, prec, nrefine, conjugate);
    func_a.CompD[0] = nullptr;
}
template <int D> void multiply(CompFunction<D> &out, FunctionTree<D, ComplexDouble> &inp_a, RepresentableFunction<D, ComplexDouble> &f, double prec, int nrefine, bool conjugate) {
    CompFunction<D> func_a(1);
    func_a.func_ptr->isreal = 0;
    func_a.func_ptr->iscomplex = 1;
    func_a.CompC[0] = &inp_a;
    multiply(out, func_a, f, prec, nrefine, conjugate);
    func_a.CompC[0] = nullptr;
}

/** @brief Compute <bra|ket> = int bra^\dag(r) * ket(r) dr.
 *
 *  Sum of component dots.
 *  Notice that the <bra| position is complex conjugated in the tree multiplication.
 *
 */
template <int D> ComplexDouble dot(CompFunction<D> bra, CompFunction<D> ket) {
    if (bra.func_ptr->conj or ket.func_ptr->conj) MSG_ABORT("Not implemented");
    ComplexDouble dotprodtot = 0.0;
    for (int comp = 0; comp < bra.Ncomp(); comp++) {
        ComplexDouble dotprod = 0.0;
        if (bra.func_ptr->data.n1[0] != ket.func_ptr->data.n1[0] and bra.func_ptr->data.n1[0] != 0 and ket.func_ptr->data.n1[0] != 0) continue;
        if (bra.isreal() and ket.isreal()) {
            dotprod += mrcpp::dot(*bra.CompD[comp], *ket.CompD[comp]);
        } else if (bra.isreal() and ket.iscomplex()) {
            dotprod += mrcpp::dot(*bra.CompD[comp], *ket.CompC[comp]);
        } else if (bra.iscomplex() and ket.isreal()) {
            dotprod += mrcpp::dot(*bra.CompC[comp], *ket.CompD[comp]);
        } else {
            dotprod += mrcpp::dot(*bra.CompC[comp], *ket.CompC[comp]);
        }
        dotprod *= std::conj(bra.func_ptr->data.c1[comp]) * ket.func_ptr->data.c1[comp];
        dotprodtot += dotprod;
    }
    return dotprodtot;
}

/** @brief Compute  <bra|ket> = int |bra^\dag(r)| * |ket(r)| dr.
 *
 *  sum of components
 */
template <int D> double node_norm_dot(CompFunction<D> bra, CompFunction<D> ket) {
    double dotprodtot = 0.0;
    for (int comp = 0; comp < bra.Ncomp(); comp++) {
        double dotprod = 0.0;
        if (bra.isreal() and ket.isreal()) {
            dotprod += mrcpp::node_norm_dot(*bra.CompD[comp], *ket.CompD[comp]);
        } else if (bra.isreal() and ket.iscomplex()) {
            MSG_ABORT("Not implemented");
        } else if (bra.iscomplex() and ket.isreal()) {
            MSG_ABORT("Not implemented");
        } else {
            dotprod += mrcpp::node_norm_dot(*bra.CompC[comp], *ket.CompC[comp]);
        }
        dotprod *= std::norm(bra.func_ptr->data.c1[comp]) * std::norm(ket.func_ptr->data.c1[comp]); // for fully complex values this does not really give the norm
        dotprodtot += dotprod;
    }
    return dotprodtot;
}

void project(CompFunction<3> &out, std::function<double(const Coord<3> &r)> f, double prec) {
    bool need_to_project = not(out.isShared()) or mpi::share_master();
    out.func_ptr->isreal = 1;
    out.func_ptr->iscomplex = 0;
    if (out.Ncomp() < 1) out.alloc(1);
    if (need_to_project) mrcpp::project<3>(prec, *out.CompD[0], f);
    mpi::share_function(out, 0, 123123, mpi::comm_share);
}

// template <int D, typename T>
void project(CompFunction<3> &out, std::function<ComplexDouble(const Coord<3> &r)> f, double prec) {
    bool need_to_project = not(out.isShared()) or mpi::share_master();
    out.func_ptr->isreal = 0;
    out.func_ptr->iscomplex = 1;
    if (out.Ncomp() < 1) out.alloc(1);
    if (need_to_project) mrcpp::project<3>(prec, *out.CompC[0], f);
    mpi::share_function(out, 0, 123123, mpi::comm_share);
}

template <int D> void project(CompFunction<D> &out, RepresentableFunction<D, double> &f, double prec) {
    bool need_to_project = not(out.isShared()) or mpi::share_master();
    out.func_ptr->isreal = 1;
    out.func_ptr->iscomplex = 0;
    if (out.Ncomp() < 1) out.alloc(1);
    if (need_to_project) build_grid(*out.CompD[0], f);
    if (need_to_project) mrcpp::project<D, double>(prec, *out.CompD[0], f);
    mpi::share_function(out, 0, 132231, mpi::comm_share);
}
template <int D> void project(CompFunction<D> &out, RepresentableFunction<D, ComplexDouble> &f, double prec) {
    bool need_to_project = not(out.isShared()) or mpi::share_master();
    out.func_ptr->isreal = 0;
    out.func_ptr->iscomplex = 1;
    if (out.Ncomp() < 1) out.alloc(1);
    if (need_to_project) build_grid(*out.CompC[0], f);
    if (need_to_project) mrcpp::project<D, ComplexDouble>(prec, *out.CompC[0], f);
    mpi::share_function(out, 0, 132231, mpi::comm_share);
}

// CompFunctionVector

CompFunctionVector::CompFunctionVector(int N)
        : std::vector<CompFunction<3>>(N) {
    for (int i = 0; i < N; i++) (*this)[i].func_ptr->rank = i;
    vecMRA = defaultCompMRA<3>;
}
void CompFunctionVector::distribute() {
    for (int i = 0; i < this->size(); i++) (*this)[i].func_ptr->rank = i;
}

/** @brief Make a linear combination of functions
 *
 * Uses "local" representation: treats one node at a time.
 * For each node, all functions are transformed simultaneously
 * by a dense matrix multiplication.
 * Phi input functions, Psi output functions
 * Phi and Psi are complex.
 */
void rotate_cplx(CompFunctionVector &Phi, const ComplexMatrix &U, CompFunctionVector &Psi, double prec) {

    // The principle of this routine is that nodes for all orbitals are rotated one by one using matrix multiplication.
    // The routine does avoid when possible to move data, but uses pointers and indices manipulation.
    // MPI version does not use OMP yet, Serial version uses OMP
    // size of input is N, size of output is M
    bool serial = mpi::wrk_size == 1; // flag for serial/MPI switch
    int N = Phi.size();
    int M = Psi.size();
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < 4; j++) delete Psi[i].CompD[j];
        Psi[i].func_ptr->isreal = 0;
        Psi[i].func_ptr->iscomplex = 1;
    }
    for (int i = 0; i < N; i++) {
        if (Phi[i].func_ptr->conj) MSG_ABORT("Conjugaison not implemneted for rotations");
    }
    if (U.rows() < N) MSG_ABORT("Incompatible number of rows for U matrix");
    if (U.cols() < M) MSG_ABORT("Incompatible number of columns for U matrix");

    // 1) make union tree without coefficients. Note that the ref tree is always real (in fact it has no coeff)
    FunctionTree<3> refTree(*Phi.vecMRA);
    mpi::allreduce_Tree_noCoeff(refTree, Phi, mpi::comm_wrk);

    int sizecoeff = (1 << refTree.getDim()) * refTree.getKp1_d();
    int sizecoeffW = ((1 << refTree.getDim()) - 1) * refTree.getKp1_d();
    std::vector<double> scalefac_ref;
    std::vector<double *> coeffVec_ref; // not used!
    std::vector<int> indexVec_ref;      // serialIx of the nodes
    std::vector<int> parindexVec_ref;   // serialIx of the parent nodes
    int max_ix;
    // get a list of all nodes in union tree, identified by their serialIx indices
    refTree.makeCoeffVector(coeffVec_ref, indexVec_ref, parindexVec_ref, scalefac_ref, max_ix, refTree);
    int max_n = indexVec_ref.size();

    for (int j = 0; j < N; j++) {
        if (!mpi::my_func(j)) continue;
        if (Phi[j].isreal()) MSG_ABORT("This function only use complex input");
    }

    for (int i = 0; i < M; i++) {
        Psi[i].func_ptr->data.isreal = 0;
        Psi[i].func_ptr->data.iscomplex = 1;
    }

    // 3) In the serial case we store the coeff pointers in coeffVec. In the mpi case the coeff are stored in the bank

    BankAccount nodesPhi;     // to put the original nodes
    BankAccount nodesRotated; // to put the rotated nodes

    // used for serial only:
    std::vector<std::vector<ComplexDouble *>> coeffVec(N);
    std::vector<std::vector<int>> indexVec(N);   // serialIx of the nodes
    std::map<int, std::vector<int>> node2orbVec; // for each node index, gives a vector with the indices of the orbitals using this node
    std::vector<std::map<int, int>> orb2node(N); // for a given orbital and a given node, gives the node index in the
                                                 // orbital given the node index in the reference tree
    if (serial) {
        // make list of all coefficients (coeffVec), and their reference indices (indexVec)
        std::vector<int> parindexVec; // serialIx of the parent nodes
        std::vector<double> scalefac;
        for (int j = 0; j < N; j++) {
            // make vector with all coef pointers and their indices in the union grid
            Phi[j].complex().makeCoeffVector(coeffVec[j], indexVec[j], parindexVec, scalefac, max_ix, refTree);
            // make a map that gives j from indexVec
            int orb_node_ix = 0;
            for (int ix : indexVec[j]) {
                orb2node[j][ix] = orb_node_ix++;
                if (ix < 0) continue;
                node2orbVec[ix].push_back(j);
            }
        }
    } else { // MPI case
        // send own nodes to bank, identifying them through the serialIx of refTree
        save_nodes(Phi, refTree, nodesPhi);
        mpi::barrier(mpi::comm_wrk); // required for now, as the blockdata functionality has no queue yet.
    }

    // 4) rotate all the nodes
    IntMatrix split_serial;                                 // in the serial case all split are stored in one array
    std::vector<std::vector<ComplexDouble *>> coeffpVec(M); // to put pointers to the rotated coefficient for each orbital in serial case
    std::vector<std::map<int, int>> ix2coef(M);             // to find the index in for example rotCoeffVec[] corresponding to a serialIx
    int csize;                                              // size of the current coefficients (different for roots and branches)
    std::vector<ComplexMatrix> rotatedCoeffVec;             // just to ensure that the data from rotatedCoeff is not deleted, since we point to it.
                                                            // j indices are for unrotated orbitals, i indices are for rotated orbitals
    if (serial) {
        std::map<int, int> ix2coef_ref; // to find the index n corresponding to a serialIx
        split_serial.resize(M, max_n);  // not use in the MPI case
        for (int n = 0; n < max_n; n++) {
            int node_ix = indexVec_ref[n]; // SerialIx for this node in the reference tree
            ix2coef_ref[node_ix] = n;
            for (int i = 0; i < M; i++) split_serial(i, n) = 1;
        }
        std::vector<int> nodeReady(max_n, 0); // To indicate to OMP threads that the parent is ready (for splits)
                                              // assumes the nodes are ordered such that parent are treated before children. BFS or DFS ok.
                                              // NB: the n must be traversed approximately in right order: Thread n may have to wait until som other preceding
                                              // n is finished.
#pragma omp parallel for schedule(dynamic)
        for (int n = 0; n < max_n; n++) {
            int csize;
            int node_ix = indexVec_ref[n]; // SerialIx for this node in the reference tree
            // 4a) make a dense contiguous matrix with the coefficient from all the orbitals using node n
            std::vector<int> orbjVec; // to remember which orbital correspond to each orbVec.size();
            if (node2orbVec[node_ix].size() <= 0) continue;
            csize = sizecoeffW;
            if (parindexVec_ref[n] < 0) csize = sizecoeff; // for root nodes we include scaling coeff

            int shift = sizecoeff - sizecoeffW; // to copy only wavelet part
            if (parindexVec_ref[n] < 0) shift = 0;
            ComplexMatrix coeffBlock(csize, node2orbVec[node_ix].size());
            for (int j : node2orbVec[node_ix]) { // loop over indices of the orbitals using this node
                int orb_node_ix = orb2node[j][node_ix];
                for (int k = 0; k < csize; k++) coeffBlock(k, orbjVec.size()) = coeffVec[j][orb_node_ix][k + shift];
                orbjVec.push_back(j);
            }

            // 4b) make a list of rotated orbitals needed for this node
            // OMP must wait until parent is ready
            while (parindexVec_ref[n] >= 0 and nodeReady[ix2coef_ref[parindexVec_ref[n]]] == 0) {
#pragma omp flush
            };

            std::vector<int> orbiVec;
            for (int i = 0; i < M; i++) {                                                                        // loop over all rotated orbitals
                if (parindexVec_ref[n] >= 0 and split_serial(i, ix2coef_ref[parindexVec_ref[n]]) == 0) continue; // parent node has too small wavelets
                orbiVec.push_back(i);
            }

            // 4c) rotate this node
            ComplexMatrix Un(orbjVec.size(), orbiVec.size()); // chunk of U, with reorganized indices
            for (int i = 0; i < orbiVec.size(); i++) {        // loop over rotated orbitals
                for (int j = 0; j < orbjVec.size(); j++) { Un(j, i) = U(orbjVec[j], orbiVec[i]); }
            }
            ComplexMatrix rotatedCoeff(csize, orbiVec.size());
            // HERE IT HAPPENS!
            // TODO: conjugaison
            rotatedCoeff.noalias() = coeffBlock * Un; // Matrix mutiplication

            // 4d) store and make rotated node pointers
            // for now we allocate in buffer, in future could be directly allocated in the final trees
            double thres = prec * prec * scalefac_ref[n] * scalefac_ref[n];
            // make all norms:
            for (int i = 0; i < orbiVec.size(); i++) {
                // check if parent must be split
                if (parindexVec_ref[n] == -1 or split_serial(orbiVec[i], ix2coef_ref[parindexVec_ref[n]])) {
                    // mark this node for this orbital for later split
#pragma omp critical
                    {
                        ix2coef[orbiVec[i]][node_ix] = coeffpVec[orbiVec[i]].size();
                        coeffpVec[orbiVec[i]].push_back(&(rotatedCoeff(0, i))); // list of coefficient pointers
                    }
                    // check norms for split
                    double wnorm = 0.0; // rotatedCoeff(k, i) is already in cache here
                    int kstart = 0;
                    if (parindexVec_ref[n] < 0) kstart = sizecoeff - sizecoeffW; // do not include scaling, even for roots
                    for (int k = kstart; k < csize; k++) wnorm += std::real(rotatedCoeff(k, i) * std::conj(rotatedCoeff(k, i)));
                    if (thres < wnorm or prec < 0)
                        split_serial(orbiVec[i], n) = 1;
                    else
                        split_serial(orbiVec[i], n) = 0;
                } else {
                    ix2coef[orbiVec[i]][node_ix] = max_n + 1; // should not be used
                    split_serial(orbiVec[i], n) = 0;          // do not split if parent does not need to be split
                }
            }
            nodeReady[n] = 1;
#pragma omp critical
            {
                // this ensures that rotatedCoeff is not deleted, when getting out of scope
                rotatedCoeffVec.push_back(std::move(rotatedCoeff));
            }
        }
    } else { // MPI case

        // TODO? rotate in bank, so that we do not get and put. Requires clever handling of splits.
        std::vector<double> split(M, -1.0);    // which orbitals need splitting (at a given node). For now double for compatibilty with bank
        std::vector<double> needsplit(M, 1.0); // which orbitals need splitting
        BankAccount nodeSplits;
        mpi::barrier(mpi::comm_wrk); // required for now, as the blockdata functionality has no queue yet.

        ComplexMatrix coeffBlock(sizecoeff, N);
        max_ix++; // largest node index + 1. to store rotated orbitals with different id
        TaskManager tasks(max_n);
        for (int nn = 0; nn < max_n; nn++) {
            int n = tasks.next_task();
            if (n < 0) break;
            double thres = prec * prec * scalefac_ref[n] * scalefac_ref[n];
            // 4a) make list of orbitals that should split the parent node, i.e. include this node
            int parentid = parindexVec_ref[n];
            if (parentid == -1) {
                // root node, split if output needed
                for (int i = 0; i < M; i++) { split[i] = 1.0; }
                csize = sizecoeff;
            } else {
                // note that it will wait until data is available
                nodeSplits.get_data(parentid, M, split.data());
                csize = sizecoeffW;
            }
            std::vector<int> orbiVec;
            std::vector<int> orbjVec;
            for (int i = 0; i < M; i++) {     // loop over rotated orbitals
                if (split[i] < 0.0) continue; // parent node has too small wavelets
                orbiVec.push_back(i);
            }

            // 4b) rotate this node
            ComplexMatrix coeffBlock(csize, N); // largest possible used size
            nodesPhi.get_nodeblock(indexVec_ref[n], coeffBlock.data(), orbjVec);
            coeffBlock.conservativeResize(Eigen::NoChange, orbjVec.size()); // keep only used part

            // chunk of U, with reorganized indices and separate blocks for real and imag:
            ComplexMatrix Un(orbjVec.size(), orbiVec.size());
            ComplexMatrix rotatedCoeff(csize, orbiVec.size());

            for (int i = 0; i < orbiVec.size(); i++) {     // loop over included rotated real and imag part of orbitals
                for (int j = 0; j < orbjVec.size(); j++) { // loop over input orbital, possibly imaginary parts
                    Un(j, i) = U(orbjVec[j], orbiVec[i]);
                }
            }

            // HERE IT HAPPENS
            // TODO conjugaison
            rotatedCoeff.noalias() = coeffBlock * Un; // Matrix mutiplication

            // 3c) find which orbitals need to further refine this node, and store rotated node (after each other while
            // in cache).
            for (int i = 0; i < orbiVec.size(); i++) { // loop over rotated orbitals
                needsplit[orbiVec[i]] = -1.0;          // default, do not split
                // check if this node/orbital needs further refinement
                double wnorm = 0.0;
                int kwstart = csize - sizecoeffW; // do not include scaling
                for (int k = kwstart; k < csize; k++) wnorm += std::real(rotatedCoeff.col(i)[k] * std::conj(rotatedCoeff.col(i)[k]));
                if (thres < wnorm or prec < 0) needsplit[orbiVec[i]] = 1.0;
                nodesRotated.put_nodedata(orbiVec[i], indexVec_ref[n] + max_ix, csize, rotatedCoeff.col(i).data());
            }
            nodeSplits.put_data(indexVec_ref[n], M, needsplit.data());
        }
        mpi::barrier(mpi::comm_wrk); // wait until all rotated nodes are ready
    }

    // 5) reconstruct trees using rotated nodes.

    // only serial case can use OMP, because MPI cannot be used by threads
    if (serial) {
        // OMP parallelized, but does not scale well, because the total memory bandwidth is a bottleneck. (the main
        // operation is writing the coefficient into the tree)

#pragma omp parallel for schedule(static)
        for (int j = 0; j < M; j++) {
            if (coeffpVec[j].size() == 0) continue;
            Psi[j].alloc(1); // All data is stored in coeffpVec[j]
            Psi[j].complex().makeTreefromCoeff(refTree, coeffpVec[j], ix2coef[j], prec);
        }
    } else { // MPI case
        for (int j = 0; j < M; j++) {
            if (not mpi::my_func(j)) continue;
            // traverse possible nodes, and stop descending when norm is zero (leaf in out[j])
            std::vector<ComplexDouble *> coeffpVec; //
            std::map<int, int> ix2coef;             // to find the index in coeffVec[] corresponding to a serialIx
            int ix = 0;
            std::vector<ComplexDouble *> pointerstodelete; // list of temporary arrays to clean up
            for (int ibank = 0; ibank < mpi::bank_size; ibank++) {
                std::vector<int> nodeidVec;
                ComplexDouble *dataVec; // will be allocated by bank
                nodesRotated.get_orbblock(j, dataVec, nodeidVec, ibank);
                if (nodeidVec.size() > 0) pointerstodelete.push_back(dataVec);
                int shift = 0;
                for (int n = 0; n < nodeidVec.size(); n++) {
                    assert(nodeidVec[n] - max_ix >= 0);                // unrotated nodes have been deleted
                    assert(ix2coef.count(nodeidVec[n] - max_ix) == 0); // each nodeid treated once
                    ix2coef[nodeidVec[n] - max_ix] = ix++;
                    csize = sizecoeffW;
                    if (parindexVec_ref[nodeidVec[n] - max_ix] < 0) csize = sizecoeff;
                    coeffpVec.push_back(&dataVec[shift]); // list of coeff pointers
                    shift += csize;
                }
            }

            Psi[j].alloc(1);
            Psi[j].complex().makeTreefromCoeff(refTree, coeffpVec, ix2coef, prec);

            for (ComplexDouble *p : pointerstodelete) delete[] p;
            pointerstodelete.clear();
        }
    }
}

/** @brief Make a linear combination of functions
 *
 * Uses "local" representation: treats one node at a time.
 * For each node, all functions are transformed simultaneously
 * by a dense matrix multiplication.
 * Phi input functions, Psi output functions
 *
 */
void rotate(CompFunctionVector &Phi, const ComplexMatrix &U, CompFunctionVector &Psi, double prec) {

    if (Phi[0].iscomplex()) {
        rotate_cplx(Phi, U, Psi, prec);
        return;
    }

    // The principle of this routine is that nodes are rotated one by one using matrix multiplication.
    // The routine does avoid when possible to move data, but uses pointers and indices manipulation.
    // MPI version does not use OMP yet, Serial version uses OMP
    // size of input is N, size of output is M
    int N = Phi.size();
    int M = Psi.size();
    if (U.rows() < N) MSG_ABORT("Incompatible number of rows for U matrix");
    if (U.cols() < M) MSG_ABORT("Incompatible number of columns for U matrix");

    // 1) make union tree without coefficients. Note that the ref tree is always real (in fact it has no coeff)
    FunctionTree<3> refTree(*Phi.vecMRA);
    mpi::allreduce_Tree_noCoeff(refTree, Phi, mpi::comm_wrk);

    int sizecoeff = (1 << refTree.getDim()) * refTree.getKp1_d();
    int sizecoeffW = ((1 << refTree.getDim()) - 1) * refTree.getKp1_d();
    std::vector<double> scalefac_ref;
    std::vector<double *> coeffVec_ref; // not used!
    std::vector<int> indexVec_ref;      // serialIx of the nodes
    std::vector<int> parindexVec_ref;   // serialIx of the parent nodes
    int max_ix;
    // get a list of all nodes in union tree, identified by their serialIx indices
    refTree.makeCoeffVector(coeffVec_ref, indexVec_ref, parindexVec_ref, scalefac_ref, max_ix, refTree);
    int max_n = indexVec_ref.size();
    for (int i = 0; i < M; i++) {
        Psi[i].func_ptr->data.isreal = 1;
        Psi[i].func_ptr->data.iscomplex = 0;
    }

    // 3) In the serial case we store the coeff pointers in coeffVec. In the mpi case the coeff are stored in the bank

    bool serial = mpi::wrk_size == 1; // flag for serial/MPI switch
    BankAccount nodesPhi;             // to put the original nodes
    BankAccount nodesRotated;         // to put the rotated nodes

    // used for serial only:
    std::vector<std::vector<double *>> coeffVec(N);
    std::vector<std::vector<int>> indexVec(N);   // serialIx of the nodes
    std::map<int, std::vector<int>> node2orbVec; // for each node index, gives a vector with the indices of the orbitals using this node
    std::vector<std::map<int, int>> orb2node(N); // for a given orbital and a given node, gives the node index in the
                                                 // orbital given the node index in the reference tree
    if (serial) {

        // make list of all coefficients (coeffVec), and their reference indices (indexVec)
        std::vector<int> parindexVec; // serialIx of the parent nodes
        std::vector<double> scalefac;
        for (int j = 0; j < N; j++) {
            // make vector with all coef pointers and their indices in the union grid
            Phi[j].real().makeCoeffVector(coeffVec[j], indexVec[j], parindexVec, scalefac, max_ix, refTree);
            // make a map that gives j from indexVec
            int orb_node_ix = 0;
            for (int ix : indexVec[j]) {
                orb2node[j][ix] = orb_node_ix++;
                if (ix < 0) continue;
                node2orbVec[ix].push_back(j);
            }
        }
    } else { // MPI case
        // send own nodes to bank, identifying them through the serialIx of refTree
        save_nodes(Phi, refTree, nodesPhi);
        mpi::barrier(mpi::comm_wrk); // required for now, as the blockdata functionality has no queue yet.
    }

    // 4) rotate all the nodes
    IntMatrix split_serial;                          // in the serial case all split are stored in one array
    std::vector<std::vector<double *>> coeffpVec(M); // to put pointers to the rotated coefficient for each orbital in serial case
    std::vector<std::map<int, int>> ix2coef(M);      // to find the index in for example rotCoeffVec[] corresponding to a serialIx
    int csize;                                       // size of the current coefficients (different for roots and branches)
    std::vector<DoubleMatrix> rotatedCoeffVec;       // just to ensure that the data from rotatedCoeff is not deleted, since we point to it.
    // j indices are for unrotated orbitals, i indices are for rotated orbitals
    if (serial) {
        std::map<int, int> ix2coef_ref; // to find the index n corresponding to a serialIx
        split_serial.resize(M, max_n);  // not use in the MPI case
        for (int n = 0; n < max_n; n++) {
            int node_ix = indexVec_ref[n]; // SerialIx for this node in the reference tree
            ix2coef_ref[node_ix] = n;
            for (int i = 0; i < M; i++) split_serial(i, n) = 1;
        }

        std::vector<int> nodeReady(max_n, 0); // To indicate to OMP threads that the parent is ready (for splits)

        // assumes the nodes are ordered such that parent are treated before children. BFS or DFS ok.
        // NB: the n must be traversed approximately in right order: Thread n may have to wait until som other preceding
        // n is finished.
#pragma omp parallel for schedule(dynamic)
        for (int n = 0; n < max_n; n++) {
            int csize;
            int node_ix = indexVec_ref[n]; // SerialIx for this node in the reference tree
            // 4a) make a dense contiguous matrix with the coefficient from all the orbitals using node n
            std::vector<int> orbjVec; // to remember which orbital correspond to each orbVec.size();
            if (node2orbVec[node_ix].size() <= 0) continue;
            csize = sizecoeffW;
            if (parindexVec_ref[n] < 0) csize = sizecoeff; // for root nodes we include scaling coeff

            int shift = sizecoeff - sizecoeffW; // to copy only wavelet part
            if (parindexVec_ref[n] < 0) shift = 0;
            DoubleMatrix coeffBlock(csize, node2orbVec[node_ix].size());
            for (int j : node2orbVec[node_ix]) { // loop over indices of the orbitals using this node
                int orb_node_ix = orb2node[j][node_ix];
                for (int k = 0; k < csize; k++) coeffBlock(k, orbjVec.size()) = coeffVec[j][orb_node_ix][k + shift];
                orbjVec.push_back(j);
            }

            // 4b) make a list of rotated orbitals needed for this node
            // OMP must wait until parent is ready
            while (parindexVec_ref[n] >= 0 and nodeReady[ix2coef_ref[parindexVec_ref[n]]] == 0) {
#pragma omp flush
            };

            std::vector<int> orbiVec;
            for (int i = 0; i < M; i++) {                                                                        // loop over all rotated orbitals
                if (parindexVec_ref[n] >= 0 and split_serial(i, ix2coef_ref[parindexVec_ref[n]]) == 0) continue; // parent node has too small wavelets
                orbiVec.push_back(i);
            }

            // 4c) rotate this node
            DoubleMatrix Un(orbjVec.size(), orbiVec.size()); // chunk of U, with reorganized indices
            for (int i = 0; i < orbiVec.size(); i++) {       // loop over rotated orbitals
                for (int j = 0; j < orbjVec.size(); j++) { Un(j, i) = std::real(U(orbjVec[j], orbiVec[i])); }
            }
            DoubleMatrix rotatedCoeff(csize, orbiVec.size());
            // HERE IT HAPPENS!
            rotatedCoeff.noalias() = coeffBlock * Un; // Matrix mutiplication

            // 4d) store and make rotated node pointers
            // for now we allocate in buffer, in future could be directly allocated in the final trees
            double thres = prec * prec * scalefac_ref[n] * scalefac_ref[n];
            // make all norms:
            for (int i = 0; i < orbiVec.size(); i++) {
                // check if parent must be split
                if (parindexVec_ref[n] == -1 or split_serial(orbiVec[i], ix2coef_ref[parindexVec_ref[n]])) {
                    // mark this node for this orbital for later split
#pragma omp critical
                    {
                        ix2coef[orbiVec[i]][node_ix] = coeffpVec[orbiVec[i]].size();
                        coeffpVec[orbiVec[i]].push_back(&(rotatedCoeff(0, i))); // list of coefficient pointers
                    }
                    // check norms for split
                    double wnorm = 0.0; // rotatedCoeff(k, i) is already in cache here
                    int kstart = 0;
                    if (parindexVec_ref[n] < 0) kstart = sizecoeff - sizecoeffW; // do not include scaling, even for roots
                    for (int k = kstart; k < csize; k++) wnorm += rotatedCoeff(k, i) * rotatedCoeff(k, i);
                    if (thres < wnorm or prec < 0)
                        split_serial(orbiVec[i], n) = 1;
                    else
                        split_serial(orbiVec[i], n) = 0;
                } else {
                    ix2coef[orbiVec[i]][node_ix] = max_n + 1; // should not be used
                    split_serial(orbiVec[i], n) = 0;          // do not split if parent does not need to be split
                }
            }
            nodeReady[n] = 1;
#pragma omp critical
            {
                // this ensures that rotatedCoeff is not deleted, when getting out of scope
                rotatedCoeffVec.push_back(std::move(rotatedCoeff));
            }
        }
    } else { // MPI case

        // TODO? rotate in bank, so that we do not get and put. Requires clever handling of splits.
        std::vector<double> split(M, -1.0);    // which orbitals need splitting (at a given node). For now double for compatibilty with bank
        std::vector<double> needsplit(M, 1.0); // which orbitals need splitting
        BankAccount nodeSplits;
        mpi::barrier(mpi::comm_wrk); // required for now, as the blockdata functionality has no queue yet.

        DoubleMatrix coeffBlock(sizecoeff, N);
        max_ix++; // largest node index + 1. to store rotated orbitals with different id
        TaskManager tasks(max_n);
        for (int nn = 0; nn < max_n; nn++) {
            int n = tasks.next_task();
            if (n < 0) break;
            double thres = prec * prec * scalefac_ref[n] * scalefac_ref[n];
            // 4a) make list of orbitals that should split the parent node, i.e. include this node
            int parentid = parindexVec_ref[n];
            if (parentid == -1) {
                // root node, split if output needed
                for (int i = 0; i < M; i++) { split[i] = 1.0; }
                csize = sizecoeff;
            } else {
                // note that it will wait until data is available
                nodeSplits.get_data(parentid, M, split.data());
                csize = sizecoeffW;
            }
            std::vector<int> orbiVec;
            std::vector<int> orbjVec;
            for (int i = 0; i < M; i++) {     // loop over rotated orbitals
                if (split[i] < 0.0) continue; // parent node has too small wavelets
                orbiVec.push_back(i);
            }

            // 4b) rotate this node
            DoubleMatrix coeffBlock(csize, N); // largest possible used size
            nodesPhi.get_nodeblock(indexVec_ref[n], coeffBlock.data(), orbjVec);
            coeffBlock.conservativeResize(Eigen::NoChange, orbjVec.size()); // keep only used part

            // chunk of U, with reorganized indices and separate blocks for real and imag:
            DoubleMatrix Un(orbjVec.size(), orbiVec.size());
            DoubleMatrix rotatedCoeff(csize, orbiVec.size());

            for (int i = 0; i < orbiVec.size(); i++) {     // loop over included rotated real and imag part of orbitals
                for (int j = 0; j < orbjVec.size(); j++) { // loop over input orbital, possibly imaginary parts
                    Un(j, i) = std::real(U(orbjVec[j], orbiVec[i]));
                }
            }

            // HERE IT HAPPENS
            rotatedCoeff.noalias() = coeffBlock * Un; // Matrix mutiplication

            // 3c) find which orbitals need to further refine this node, and store rotated node (after each other while
            // in cache).
            for (int i = 0; i < orbiVec.size(); i++) { // loop over rotated orbitals
                needsplit[orbiVec[i]] = -1.0;          // default, do not split
                // check if this node/orbital needs further refinement
                double wnorm = 0.0;
                int kwstart = csize - sizecoeffW; // do not include scaling
                for (int k = kwstart; k < csize; k++) wnorm += rotatedCoeff.col(i)[k] * rotatedCoeff.col(i)[k];
                if (thres < wnorm or prec < 0) needsplit[orbiVec[i]] = 1.0;
                nodesRotated.put_nodedata(orbiVec[i], indexVec_ref[n] + max_ix, csize, rotatedCoeff.col(i).data());
            }
            nodeSplits.put_data(indexVec_ref[n], M, needsplit.data());
        }
        mpi::barrier(mpi::comm_wrk); // wait until all rotated nodes are ready
    }

    // 5) reconstruct trees using rotated nodes.

    // only serial case can use OMP, because MPI cannot be used by threads
    if (serial) {
        // OMP parallelized, but does not scale well, because the total memory bandwidth is a bottleneck. (the main
        // operation is writing the coefficient into the tree)

#pragma omp parallel for schedule(static)
        for (int j = 0; j < M; j++) {
            if (coeffpVec[j].size() == 0) continue;
            Psi[j].alloc(1);
            Psi[j].real().clear();
            Psi[j].real().makeTreefromCoeff(refTree, coeffpVec[j], ix2coef[j], prec);
        }

    } else { // MPI case

        for (int j = 0; j < M; j++) {
            if (not mpi::my_func(j)) continue;
            // traverse possible nodes, and stop descending when norm is zero (leaf in out[j])
            std::vector<double *> coeffpVec; //
            std::map<int, int> ix2coef;      // to find the index in coeffVec[] corresponding to a serialIx
            int ix = 0;
            std::vector<double *> pointerstodelete; // list of temporary arrays to clean up
            for (int ibank = 0; ibank < mpi::bank_size; ibank++) {
                std::vector<int> nodeidVec;
                double *dataVec; // will be allocated by bank
                nodesRotated.get_orbblock(j, dataVec, nodeidVec, ibank);
                if (nodeidVec.size() > 0) pointerstodelete.push_back(dataVec);
                int shift = 0;
                for (int n = 0; n < nodeidVec.size(); n++) {
                    assert(nodeidVec[n] - max_ix >= 0);                // unrotated nodes have been deleted
                    assert(ix2coef.count(nodeidVec[n] - max_ix) == 0); // each nodeid treated once
                    ix2coef[nodeidVec[n] - max_ix] = ix++;
                    csize = sizecoeffW;
                    if (parindexVec_ref[nodeidVec[n] - max_ix] < 0) csize = sizecoeff;
                    coeffpVec.push_back(&dataVec[shift]); // list of coeff pointers
                    shift += csize;
                }
            }
            Psi[j].alloc(1);
            Psi[j].real().makeTreefromCoeff(refTree, coeffpVec, ix2coef, prec);

            for (double *p : pointerstodelete) delete[] p;
            pointerstodelete.clear();
        }
    }
}

void rotate(CompFunctionVector &Phi, const ComplexMatrix &U, double prec) {
    rotate(Phi, U, Phi, prec);
    return;
}

/** @brief Save all nodes in bank; identify them using serialIx from refTree
 * shift is a shift applied in the id
 */
void save_nodes(CompFunctionVector &Phi, FunctionTree<3> &refTree, BankAccount &account, int sizes) {
    int sizecoeff = (1 << refTree.getDim()) * refTree.getKp1_d();
    int sizecoeffW = ((1 << refTree.getDim()) - 1) * refTree.getKp1_d();
    int max_nNodes = refTree.getNNodes();
    std::vector<double *> coeffVec;
    std::vector<ComplexDouble *> coeffVec_cplx;
    std::vector<double> scalefac;
    std::vector<int> indexVec;    // SerialIx of the node in refOrb
    std::vector<int> parindexVec; // SerialIx of the parent node
    int N = Phi.size();
    int max_ix;
    for (int j = 0; j < N; j++) {
        if (not mpi::my_func(j)) continue;
        // make vector with all coef address and their index in the union grid
        if (Phi[j].isreal()) {
            Phi[j].real().makeCoeffVector(coeffVec, indexVec, parindexVec, scalefac, max_ix, refTree);
            int max_n = indexVec.size();
            // send node coefs from Phi[j] to bank
            // except for the root nodes, only wavelets are sent
            for (int i = 0; i < max_n; i++) {
                if (indexVec[i] < 0) continue; // nodes that are not in refOrb
                int csize = sizecoeffW;
                if (parindexVec[i] < 0) csize = sizecoeff;
                if (sizes > 0) { // fixed size
                    account.put_nodedata(j, indexVec[i], sizes, coeffVec[i]);
                } else {
                    account.put_nodedata(j, indexVec[i], csize, &(coeffVec[i][sizecoeff - csize]));
                }
            }
        }
        // Complex components
        if (Phi[j].iscomplex()) {
            Phi[j].complex().makeCoeffVector(coeffVec_cplx, indexVec, parindexVec, scalefac, max_ix, refTree);
            int max_n = indexVec.size();
            // send node coefs from Phi[j] to bank
            for (int i = 0; i < max_n; i++) {
                if (indexVec[i] < 0) continue; // nodes that are not in refOrb
                // NB: the identifier (indexVec[i]) must be shifted for not colliding with the nodes from the real part
                int csize = sizecoeffW;
                if (parindexVec[i] < 0) csize = sizecoeff;
                if (sizes > 0) { // fixed size
                    account.put_nodedata(j, indexVec[i], sizes, coeffVec_cplx[i]);
                } else {
                    account.put_nodedata(j, indexVec[i], csize, &(coeffVec_cplx[i][sizecoeff - csize]));
                }
            }
        }
    }
}

/** @brief Multiply all orbitals with a function
 *
 * @param Phi: orbitals to multiply
 * @param f  : function to multiply
 *
 * Computes the product of each orbital with a function
 * in parallel using a local representation.
 * Input trees are extended by one scale at most.
 */
CompFunctionVector multiply(CompFunctionVector &Phi, RepresentableFunction<3> &f, double prec, CompFunction<3> *Func, int nrefine, bool all) {
    int N = Phi.size();
    const int D = 3;
    bool serial = mpi::wrk_size == 1; // flag for serial/MPI switch
    // 1a) extend grid where f is large (around nuclei)
    // TODO: do it in save_nodes + refTree, only saving the extra nodes, without keeping them permanently. Or refine refTree?

    for (int i = 0; i < N; i++) {
        if (!mpi::my_func(i)) continue;
        int irefine = 0;
        while (Phi[i].isreal() and irefine < nrefine and refine_grid(Phi[i].real(), f) > 0) irefine++;
        if (Phi[i].iscomplex()) MSG_ABORT("Not yet implemented");
        irefine = 0;
        //        while (Phi[i].iscomplex() and irefine < nrefine and refine_grid(Phi[i].complex(), f) > 0) irefine++;
    }

    // 1b) make union tree without coefficients
    FunctionTree<D> refTree(*Phi.vecMRA);
    // refine_grid(refTree, f); //to test
    mpi::allreduce_Tree_noCoeff(refTree, Phi, mpi::comm_wrk);

    int kp1 = refTree.getKp1();
    int kp1_d = refTree.getKp1_d();
    int nCoefs = refTree.getTDim() * kp1_d;

    IntVector PsihasReIm = IntVector::Zero(2);
    for (int i = 0; i < N; i++) {
        if (!mpi::my_func(i)) continue;
        PsihasReIm[0] = (Phi[i].hasReal()) ? 1 : 0;
        PsihasReIm[1] = (Phi[i].hasImag()) ? 1 : 0;
    }
    mpi::allreduce_vector(PsihasReIm, mpi::comm_wrk);
    CompFunctionVector out(N);
    for (int i = 0; i < N; i++) { out[0] = Phi[i].paramCopy(); }
    if (not PsihasReIm[0] and not PsihasReIm[1]) {
        return out; // do nothing
    }

    std::vector<double> scalefac_ref;
    std::vector<double *> coeffVec_ref; // not used!
    std::vector<int> indexVec_ref;      // serialIx of the nodes
    std::vector<int> parindexVec_ref;   // serialIx of the parent nodes
    std::vector<MWNode<D> *> refNodes;  // pointers to nodes
    int max_ix;
    // get a list of all nodes in union tree, identified by their serialIx indices
    refTree.makeCoeffVector(coeffVec_ref, indexVec_ref, parindexVec_ref, scalefac_ref, max_ix, refTree, &refNodes);
    int max_n = indexVec_ref.size();
    std::map<int, int> ix2n; // for a given serialIx, give index in vectors
    for (int nn = 0; nn < max_n; nn++) ix2n[indexVec_ref[nn]] = nn;

    // 2a) send own nodes to bank, identifying them through the serialIx of refTree
    BankAccount nodesPhi;        // to put the original nodes
    BankAccount nodesMultiplied; // to put the multiplied nodes

    // used for serial only:
    std::vector<std::vector<double *>> coeffVec(N);
    std::vector<std::vector<int>> indexVec(N);   // serialIx of the nodes
    std::map<int, std::vector<int>> node2orbVec; // for each node index, gives a vector with the indices of the orbitals using this node
    std::vector<std::map<int, int>> orb2node(N); // for a given orbital and a given node, gives the node index in the
                                                 // orbital given the node index in the reference tree
    if (serial) {
        // make list of all coefficients (coeffVec), and their reference indices (indexVec)
        std::vector<int> parindexVec; // serialIx of the parent nodes
        std::vector<double> scalefac;
        for (int j = 0; j < N; j++) {
            // make vector with all coef pointers and their indices in the union grid
            if (Phi[j].hasReal()) {
                Phi[j].real().makeCoeffVector(coeffVec[j], indexVec[j], parindexVec, scalefac, max_ix, refTree);
                // make a map that gives j from indexVec
                int orb_node_ix = 0;
                for (int ix : indexVec[j]) {
                    orb2node[j][ix] = orb_node_ix++;
                    if (ix < 0) continue;
                    node2orbVec[ix].push_back(j);
                }
            }
            if (Phi[j].hasImag()) {
                Phi[j].imag().makeCoeffVector(coeffVec[j + N], indexVec[j + N], parindexVec, scalefac, max_ix, refTree);
                // make a map that gives j from indexVec
                int orb_node_ix = 0;
                for (int ix : indexVec[j + N]) {
                    orb2node[j + N][ix] = orb_node_ix++;
                    if (ix < 0) continue;
                    node2orbVec[ix].push_back(j + N);
                }
            }
        }
    } else {
        save_nodes(Phi, refTree, nodesPhi, nCoefs);
        mpi::barrier(mpi::comm_wrk); // required for now, as the blockdata functionality has no queue yet.
    }

    // 2b) save Func in bank and remove its coefficients
    if (Func != nullptr and !serial) {
        // put Func in local representation if not already done
        if (!Func->real().isLocal) { Func->real().saveNodesAndRmCoeff(); }
    }

    // 3) mutiply for each node
    std::vector<std::vector<double *>> coeffpVec(N); // to put pointers to the multiplied coefficient for each orbital in serial case
    std::vector<DoubleMatrix> multipliedCoeffVec;    // just to ensure that the data from multipliedCoeff is not deleted, since we point to it.
    std::vector<std::map<int, int>> ix2coef(N);      // to find the index in for example rotCoeffVec[] corresponding to a serialIx
    DoubleVector NODEP = DoubleVector::Zero(nCoefs);
    DoubleVector NODEF = DoubleVector::Zero(nCoefs);

    if (serial) {
#pragma omp parallel for schedule(dynamic)
        for (int n = 0; n < max_n; n++) {
            MWNode<D> node(*(refNodes[n]), false);
            int node_ix = indexVec_ref[n]; // SerialIx for this node in the reference tree

            // 3a) make values for f at this node
            // 3a1) get coordinates of quadrature points for this node
            Eigen::MatrixXd pts; // Eigen::Zero(D, nCoefs);
            double fval[nCoefs];
            Coord<D> r;
            double *originalCoef = nullptr;
            MWNode<3> *Fnode = nullptr;
            if (Func == nullptr) {
                node.getExpandedChildPts(pts); // TODO: use getPrimitiveChildPts (less cache).
                for (int j = 0; j < nCoefs; j++) {
                    for (int d = 0; d < D; d++) r[d] = pts(d, j); //*scaling_factor[d]?
                    fval[j] = f.evalf(r);
                }
            } else {
                Fnode = Func->real().findNode(node.getNodeIndex());
                if (Fnode == nullptr) {
                    node.getExpandedChildPts(pts); // TODO: use getPrimitiveChildPts (less cache).
                    for (int j = 0; j < nCoefs; j++) {
                        for (int d = 0; d < D; d++) r[d] = pts(d, j); //*scaling_factor[d]?
                        fval[j] = f.evalf(r);
                    }
                } else {
                    originalCoef = Fnode->getCoefs();
                    for (int j = 0; j < nCoefs; j++) fval[j] = originalCoef[j];
                    Fnode->attachCoefs(fval); // note that each thread has its own copy
                    Fnode->mwTransform(Reconstruction);
                    Fnode->cvTransform(Forward);
                }
            }
            DoubleMatrix multipliedCoeff(nCoefs, node2orbVec[node_ix].size());
            int i = 0;
            // 3b) fetch all orbitals at this node
            std::vector<int> orbjVec;            // to remember which orbital correspond to each orbVec.size();
            for (int j : node2orbVec[node_ix]) { // loop over indices of the orbitals using this node
                int orb_node_ix = orb2node[j][node_ix];
                orbjVec.push_back(j);
                for (int k = 0; k < nCoefs; k++) multipliedCoeff(k, i) = coeffVec[j][orb_node_ix][k];
                // 3c) transform to grid
                node.attachCoefs(&(multipliedCoeff(0, i)));
                node.mwTransform(Reconstruction);
                node.cvTransform(Forward);
                // 3d) multiply
                for (int k = 0; k < nCoefs; k++) multipliedCoeff(k, i) *= fval[k]; // replace by Matrix vector multiplication?
                // 3e) transform back to mw
                node.cvTransform(Backward);
                node.mwTransform(Compression);
                i++;
            }
            if (Func != nullptr and originalCoef != nullptr) {
                // restablish original values
                Fnode->attachCoefs(originalCoef);
            }

            // 3f) save multiplied nodes
            for (int i = 0; i < orbjVec.size(); i++) {
#pragma omp critical
                {
                    ix2coef[orbjVec[i]][node_ix] = coeffpVec[orbjVec[i]].size();
                    coeffpVec[orbjVec[i]].push_back(&(multipliedCoeff(0, i))); // list of coefficient pointers
                }
            }
#pragma omp critical
            {
                // this ensures that multipliedCoeff is not deleted, when getting out of scope
                multipliedCoeffVec.push_back(std::move(multipliedCoeff));
            }
            node.attachCoefs(nullptr); // to avoid deletion of valid multipliedCoeff by destructor
        }
    } else {
        // MPI
        int count1 = 0;
        int count2 = 0;
        TaskManager tasks(max_n);
        for (int nn = 0; nn < max_n; nn++) {
            int n = tasks.next_task();
            if (n < 0) break;
            MWNode<D> node(*(refNodes[n]), false);
            // 3a) make values for f
            // 3a1) get coordinates of quadrature points for this node
            Eigen::MatrixXd pts;           // Eigen::Zero(D, nCoefs);
            node.getExpandedChildPts(pts); // TODO: use getPrimitiveChildPts (less cache).
            double fval[nCoefs];
            Coord<D> r;
            MWNode<D> Fnode(*(refNodes[n]), false);
            if (Func == nullptr) {
                for (int j = 0; j < nCoefs; j++) {
                    for (int d = 0; d < D; d++) r[d] = pts(d, j); //*scaling_factor[d]?
                    fval[j] = f.evalf(r);
                }
            } else {
                int nIdx = Func->real().getIx(node.getNodeIndex());
                count1++;
                if (nIdx < 0) {
                    // use the function f instead of Func
                    count2++;
                    for (int j = 0; j < nCoefs; j++) {
                        for (int d = 0; d < D; d++) r[d] = pts(d, j);
                        fval[j] = f.evalf(r);
                    }
                } else {
                    Func->real().getNodeCoeff(nIdx, fval); // fetch coef from Bank
                    Fnode.attachCoefs(fval);
                    Fnode.mwTransform(Reconstruction);
                    Fnode.cvTransform(Forward);
                }
            }

            // 3b) fetch all orbitals at this node
            DoubleMatrix coeffBlock(nCoefs, N); // largest possible used size
            std::vector<int> orbjVec;
            nodesPhi.get_nodeblock(indexVec_ref[n], coeffBlock.data(), orbjVec);
            coeffBlock.conservativeResize(Eigen::NoChange, orbjVec.size()); // keep only used part
            DoubleMatrix MultipliedCoeff(nCoefs, orbjVec.size());
            // 3c) transform to grid
            for (int j = 0; j < orbjVec.size(); j++) { // TODO: transform all j at once ?
                // TODO: select only nodes that are end nodes?
                node.attachCoefs(coeffBlock.col(j).data());
                node.mwTransform(Reconstruction);
                node.cvTransform(Forward);
                // 3d) multiply
                double *coefs = node.getCoefs();
                for (int i = 0; i < nCoefs; i++) coefs[i] *= fval[i];
                // 3e) transform back to mw
                node.cvTransform(Backward);
                node.mwTransform(Compression);
                // 3f) save multiplied nodes
                nodesMultiplied.put_nodedata(orbjVec[j], indexVec_ref[n] + max_ix, nCoefs, coefs);
            }
            node.attachCoefs(nullptr);  // to avoid deletion of valid multipliedCoeff by destructor
            Fnode.attachCoefs(nullptr); // to avoid deletion of valid multipliedCoeff by destructor
        }
        mrcpp::mpi::barrier(mrcpp::mpi::comm_wrk); // wait until everything is stored before fetching!
    }

    // 5) reconstruct trees using multiplied nodes.

    // only serial case can use OMP, because MPI cannot be used by threads
    if (serial) {
        // OMP parallelized, but does not scale well, because the total memory bandwidth is a bottleneck. (the main
        // operation is writing the coefficient into the tree)

#pragma omp parallel for schedule(static)
        for (int j = 0; j < N; j++) {
            if (j < N) {
                if (Phi[j].hasReal()) {
                    out[j].alloc(1);
                    out[j].real().clear();
                    out[j].real().makeTreefromCoeff(refTree, coeffpVec[j], ix2coef[j], -1.0, "copy");
                    // 6) reconstruct trees from end nodes
                    out[j].real().mwTransform(BottomUp);
                    out[j].real().calcSquareNorm();
                }
            } else {
                if (Phi[j].hasImag()) {
                    out[j].alloc(1);
                    out[j].imag().clear();
                    out[j].imag().makeTreefromCoeff(refTree, coeffpVec[j], ix2coef[j], -1.0, "copy");
                    out[j].imag().mwTransform(BottomUp);
                    out[j].imag().calcSquareNorm();
                }
            }
        }
    } else {
        for (int j = 0; j < N; j++) {
            if (not mpi::my_func(j) and not all) continue;
            // traverse possible nodes, and stop descending when norm is zero (leaf in out[j])
            std::vector<double *> coeffpVec; //
            std::map<int, int> ix2coef;      // to find the index in coeffVec[] corresponding to a serialIx in refTree
            int ix = 0;
            std::vector<double *> pointerstodelete; // list of temporary arrays to clean up

            for (int ibank = 0; ibank < mpi::bank_size; ibank++) {
                std::vector<int> nodeidVec;
                double *dataVec; // will be allocated by bank
                nodesMultiplied.get_orbblock(j, dataVec, nodeidVec, ibank);
                if (nodeidVec.size() > 0) pointerstodelete.push_back(dataVec);
                int shift = 0;
                for (int n = 0; n < nodeidVec.size(); n++) {
                    assert(nodeidVec[n] - max_ix >= 0);                // unmultiplied nodes have been deleted
                    assert(ix2coef.count(nodeidVec[n] - max_ix) == 0); // each nodeid treated once
                    ix2coef[nodeidVec[n] - max_ix] = ix++;
                    coeffpVec.push_back(&dataVec[shift]); // list of coeff pointers
                    shift += nCoefs;
                }
            }
            if (j < N) {
                if (Phi[j].hasReal()) {
                    out[j].alloc(1);
                    out[j].real().clear();
                    out[j].real().makeTreefromCoeff(refTree, coeffpVec, ix2coef, -1.0, "copy");
                    // 6) reconstruct trees from end nodes
                    out[j].real().mwTransform(BottomUp);
                    out[j].real().calcSquareNorm();
                    out[j].real().resetEndNodeTable();
                    // out[j].real().crop(prec, 1.0, false); //bad convergence if out is cropped
                    if (nrefine > 0) Phi[j].real().crop(prec, 1.0, false); // restablishes original Phi
                }
            } else {
                if (Phi[j].hasImag()) {
                    out[j].alloc(1);
                    out[j].imag().clear();
                    out[j].imag().makeTreefromCoeff(refTree, coeffpVec, ix2coef, -1.0, "copy");
                    out[j].imag().mwTransform(BottomUp);
                    out[j].imag().calcSquareNorm();
                    // out[j].imag().crop(prec, 1.0, false);
                    if (nrefine > 0) Phi[j].imag().crop(prec, 1.0, false);
                }
            }

            for (double *p : pointerstodelete) delete[] p;
            pointerstodelete.clear();
        }
    }
    return out;
}

void SetdefaultMRA(MultiResolutionAnalysis<3> *MRA) {
    defaultCompMRA<3> = MRA;
}

ComplexVector dot(CompFunctionVector &Bra, CompFunctionVector &Ket) {
    int N = Bra.size();
    ComplexVector result = ComplexVector::Zero(N);
    for (int i = 0; i < N; i++) {
        // The bra is sent to the owner of the ket
        if (my_func(Bra[i]) != my_func(Ket[i])) { MSG_ABORT("same indices should have same ownership"); }
        result[i] = dot(Bra[i], Ket[i]);
        if (not mrcpp::mpi::my_func(i)) Bra[i].free();
    }
    mrcpp::mpi::allreduce_vector(result, mrcpp::mpi::comm_wrk);
    return result;
}

/** @brief Compute Löwdin orthonormalization matrix
 *
 * @param Phi: orbitals to orthonomalize
 *
 * Computes the inverse square root of the orbital overlap matrix S^(-1/2)
 */
ComplexMatrix calc_lowdin_matrix(CompFunctionVector &Phi) {
    ComplexMatrix S_tilde = calc_overlap_matrix(Phi);
    ComplexMatrix S_m12 = math_utils::hermitian_matrix_pow(S_tilde, -1.0 / 2.0);
    return S_m12;
}

/** @brief Orbital transformation out_j = sum_i inp_i*U_ij
 *
 * NOTE: OrbitalVector is considered a ROW vector, so rotation
 *       means matrix multiplication from the right
 *
 * MPI: Rank distribution of output vector is the same as input vector
 *
 */
ComplexMatrix calc_overlap_matrix_cplx(CompFunctionVector &BraKet) {
    int N = BraKet.size();
    ComplexMatrix S = ComplexMatrix::Zero(N, N);
    DoubleMatrix Sreal = S.real();
    MultiResolutionAnalysis<3> *mra = BraKet.vecMRA;

    // 1) make union tree without coefficients
    mrcpp::FunctionTree<3> refTree(*mra);
    mpi::allreduce_Tree_noCoeff(refTree, BraKet, mpi::comm_wrk);

    int sizecoeff = (1 << refTree.getDim()) * refTree.getKp1_d();
    int sizecoeffW = ((1 << refTree.getDim()) - 1) * refTree.getKp1_d();

    // get a list of all nodes in union grid, as defined by their indices
    std::vector<double> scalefac;
    std::vector<double *> coeffVec_ref;
    std::vector<int> indexVec_ref;    // serialIx of the nodes
    std::vector<int> parindexVec_ref; // serialIx of the parent nodes
    int max_ix;                       // largest index value (not used here)

    refTree.makeCoeffVector(coeffVec_ref, indexVec_ref, parindexVec_ref, scalefac, max_ix, refTree);
    int max_n = indexVec_ref.size();

    // only used for serial case:
    std::vector<std::vector<ComplexDouble *>> coeffVec(N);
    std::map<int, std::vector<int>> node2orbVec; // for each node index, gives a vector with the indices of the orbitals using this node
    std::vector<std::map<int, int>> orb2node(N); // for a given orbital and a given node, gives the node index in
                                                 // the orbital given the node index in the reference tree

    bool serial = mrcpp::mpi::wrk_size == 1; // flag for serial/MPI switch
    mrcpp::BankAccount nodesBraKet;

    // In the serial case we store the coeff pointers in coeffVec. In the mpi case the coeff are stored in the bank
    if (serial) {
        // 2) make list of all coefficients, and their reference indices
        // for different orbitals, indexVec will give the same index for the same node in space
        std::vector<int> parindexVec; // serialIx of the parent nodes
        std::vector<int> indexVec;    // serialIx of the nodes
        for (int j = 0; j < N; j++) {
            // make vector with all coef pointers and their indices in the union grid
            BraKet[j].complex().makeCoeffVector(coeffVec[j], indexVec, parindexVec, scalefac, max_ix, refTree);
            // make a map that gives j from indexVec
            int orb_node_ix = 0;
            for (int ix : indexVec) {
                orb2node[j][ix] = orb_node_ix++;
                if (ix < 0) continue;
                node2orbVec[ix].push_back(j);
            }
        }
    } else { // MPI case
        // 2) send own nodes to bank, identifying them through the serialIx of refTree
        save_nodes(BraKet, refTree, nodesBraKet);
        mrcpp::mpi::barrier(mrcpp::mpi::comm_wrk); // wait until everything is stored before fetching!
    }

    // 3) make dot product for all the nodes and accumulate into S
    int ibank = 0;
#pragma omp parallel if (serial)
    {
        ComplexMatrix S_omp = ComplexMatrix::Zero(N, N); // copy for each thread

#pragma omp for schedule(dynamic)
        for (int n = 0; n < max_n; n++) {
            if (n % mrcpp::mpi::wrk_size != mrcpp::mpi::wrk_rank) continue;
            int csize;
            int node_ix = indexVec_ref[n]; // SerialIx for this node in the reference tree
            std::vector<int> orbVec;       // identifies which orbitals use this node
            if (serial and node2orbVec[node_ix].size() <= 0) continue;
            if (parindexVec_ref[n] < 0)
                csize = sizecoeff;
            else
                csize = sizecoeffW;

            // In the serial case we copy the coeff coeffBlock. In the mpi case coeffBlock is provided by the bank
            if (serial) {
                int shift = sizecoeff - sizecoeffW; // to copy only wavelet part
                if (parindexVec_ref[n] < 0) shift = 0;
                ComplexMatrix coeffBlock(csize, node2orbVec[node_ix].size());
                for (int j : node2orbVec[node_ix]) { // loop over indices of the orbitals using this node
                    int orb_node_ix = orb2node[j][node_ix];
                    for (int k = 0; k < csize; k++) coeffBlock(k, orbVec.size()) = coeffVec[j][orb_node_ix][k + shift];
                    orbVec.push_back(j);
                }
                if (orbVec.size() > 0) {
                    ComplexMatrix S_temp(orbVec.size(), orbVec.size());
                    S_temp.noalias() = coeffBlock.transpose().conjugate() * coeffBlock;
                    for (int i = 0; i < orbVec.size(); i++) {
                        for (int j = 0; j < orbVec.size(); j++) {
                            if (BraKet[orbVec[i]].func_ptr->data.n1[0] != BraKet[orbVec[j]].func_ptr->data.n1[0] and BraKet[orbVec[i]].func_ptr->data.n1[0] != 0 and
                                BraKet[orbVec[j]].func_ptr->data.n1[0] != 0)
                                continue;
                            S_omp(orbVec[i], orbVec[j]) += S_temp(i, j);
                        }
                    }
                }
            } else { // MPI case
                ComplexMatrix coeffBlock(csize, N);
                nodesBraKet.get_nodeblock(indexVec_ref[n], coeffBlock.data(), orbVec);

                if (orbVec.size() > 0) {
                    ComplexMatrix S_temp(orbVec.size(), orbVec.size());
                    coeffBlock.conservativeResize(Eigen::NoChange, orbVec.size());
                    S_temp.noalias() = coeffBlock.transpose().conjugate() * coeffBlock;
                    for (int i = 0; i < orbVec.size(); i++) {
                        for (int j = 0; j < orbVec.size(); j++) {
                            if (BraKet[orbVec[i]].func_ptr->data.n1[0] != BraKet[orbVec[j]].func_ptr->data.n1[0] and BraKet[orbVec[i]].func_ptr->data.n1[0] != 0 and
                                BraKet[orbVec[j]].func_ptr->data.n1[0] != 0)
                                continue;
                            S_omp(orbVec[i], orbVec[j]) += S_temp(i, j);
                        }
                    }
                }
            }
        }
        if (serial) {
#pragma omp critical
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) { S(i, j) += S_omp(i, j); }
            }
        }
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j <= i; j++) {
            if (i != j) S(j, i) = std::conj(S(i, j)); // ensure exact symmetri
        }
    }

    // Assumes linearity: result is sum of all nodes contributions
    mrcpp::mpi::allreduce_matrix(S, mrcpp::mpi::comm_wrk);
    // multiply by CompFunction multiplicative factor

    ComplexVector Fac = ComplexVector::Zero(N);
    for (int i = 0; i < N; i++) {
        if (!mrcpp::mpi::my_func(BraKet[i])) continue;
        Fac[i] = BraKet[i].func_ptr->data.c1[0];
    }

    mrcpp::mpi::allreduce_vector(Fac, mrcpp::mpi::comm_wrk);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) { S(i, j) *= std::conj(Fac[i]) * Fac[j]; }
    }

    return S;
}
ComplexMatrix calc_overlap_matrix(CompFunctionVector &BraKet) {
    // NB: should be spinseparated at this point!
    if (BraKet[0].iscomplex()) { return calc_overlap_matrix_cplx(BraKet); }

    int N = BraKet.size();
    ComplexMatrix S = ComplexMatrix::Zero(N, N);

    MultiResolutionAnalysis<3> *mra = BraKet.vecMRA;

    // 1) make union tree without coefficients
    mrcpp::FunctionTree<3> refTree(*mra);
    mpi::allreduce_Tree_noCoeff(refTree, BraKet, mpi::comm_wrk);

    int sizecoeff = (1 << refTree.getDim()) * refTree.getKp1_d();
    int sizecoeffW = ((1 << refTree.getDim()) - 1) * refTree.getKp1_d();

    // get a list of all nodes in union grid, as defined by their indices
    std::vector<double> scalefac;
    std::vector<double *> coeffVec_ref;
    std::vector<int> indexVec_ref;    // serialIx of the nodes
    std::vector<int> parindexVec_ref; // serialIx of the parent nodes
    int max_ix;                       // largest index value (not used here)

    refTree.makeCoeffVector(coeffVec_ref, indexVec_ref, parindexVec_ref, scalefac, max_ix, refTree);
    int max_n = indexVec_ref.size();

    // only used for serial case:
    std::vector<std::vector<double *>> coeffVec(N);
    std::map<int, std::vector<int>> node2orbVec; // for each node index, gives a vector with the indices of the orbitals using this node
    std::vector<std::map<int, int>> orb2node(N); // for a given orbital and a given node, gives the node index in
                                                 // the orbital given the node index in the reference tree

    bool serial = mrcpp::mpi::wrk_size == 1; // flag for serial/MPI switch
    mrcpp::BankAccount nodesBraKet;

    // In the serial case we store the coeff pointers in coeffVec. In the mpi case the coeff are stored in the bank
    if (serial) {
        // 2) make list of all coefficients, and their reference indices
        // for different orbitals, indexVec will give the same index for the same node in space
        std::vector<int> parindexVec; // serialIx of the parent nodes
        std::vector<int> indexVec;    // serialIx of the nodes
        for (int j = 0; j < N; j++) {
            // make vector with all coef pointers and their indices in the union grid
            BraKet[j].real().makeCoeffVector(coeffVec[j], indexVec, parindexVec, scalefac, max_ix, refTree);
            // make a map that gives j from indexVec
            int orb_node_ix = 0;
            for (int ix : indexVec) {
                orb2node[j][ix] = orb_node_ix++;
                if (ix < 0) continue;
                node2orbVec[ix].push_back(j);
            }
        }
    } else { // MPI case
        // 2) send own nodes to bank, identifying them through the serialIx of refTree
        save_nodes(BraKet, refTree, nodesBraKet);
        mrcpp::mpi::barrier(mrcpp::mpi::comm_wrk); // wait until everything is stored before fetching!
    }

    // 3) make dot product for all the nodes and accumulate into S
    int ibank = 0;
#pragma omp parallel if (serial)
    {
        ComplexMatrix S_omp = ComplexMatrix::Zero(N, N); // copy for each thread

#pragma omp for schedule(dynamic)
        for (int n = 0; n < max_n; n++) {
            if (n % mrcpp::mpi::wrk_size != mrcpp::mpi::wrk_rank) continue;
            int csize;
            int node_ix = indexVec_ref[n]; // SerialIx for this node in the reference tree
            std::vector<int> orbVec;       // identifies which orbitals use this node
            if (serial and node2orbVec[node_ix].size() <= 0) continue;
            if (parindexVec_ref[n] < 0)
                csize = sizecoeff;
            else
                csize = sizecoeffW;

            // In the serial case we copy the coeff coeffBlock. In the mpi case coeffBlock is provided by the bank
            if (serial) {
                int shift = sizecoeff - sizecoeffW; // to copy only wavelet part
                if (parindexVec_ref[n] < 0) shift = 0;
                DoubleMatrix coeffBlock(csize, node2orbVec[node_ix].size());
                for (int j : node2orbVec[node_ix]) { // loop over indices of the orbitals using this node
                    int orb_node_ix = orb2node[j][node_ix];
                    for (int k = 0; k < csize; k++) coeffBlock(k, orbVec.size()) = coeffVec[j][orb_node_ix][k + shift];
                    orbVec.push_back(j);
                }
                if (orbVec.size() > 0) {
                    ComplexMatrix S_temp(orbVec.size(), orbVec.size());
                    S_temp.noalias() = coeffBlock.transpose() * coeffBlock;
                    for (int i = 0; i < orbVec.size(); i++) {
                        for (int j = 0; j < orbVec.size(); j++) {
                            if (BraKet[orbVec[i]].func_ptr->data.n1[0] != BraKet[orbVec[j]].func_ptr->data.n1[0] and BraKet[orbVec[i]].func_ptr->data.n1[0] != 0 and
                                BraKet[orbVec[j]].func_ptr->data.n1[0] != 0)
                                continue;
                            S_omp(orbVec[i], orbVec[j]) += S_temp(i, j);
                        }
                    }
                }
            } else { // MPI case
                DoubleMatrix coeffBlock(csize, N);
                nodesBraKet.get_nodeblock(indexVec_ref[n], coeffBlock.data(), orbVec);

                if (orbVec.size() > 0) {
                    DoubleMatrix S_temp(orbVec.size(), orbVec.size());
                    coeffBlock.conservativeResize(Eigen::NoChange, orbVec.size());
                    S_temp.noalias() = coeffBlock.transpose() * coeffBlock;
                    for (int i = 0; i < orbVec.size(); i++) {
                        for (int j = 0; j < orbVec.size(); j++) {
                            if (BraKet[orbVec[i]].func_ptr->data.n1[0] != BraKet[orbVec[j]].func_ptr->data.n1[0] and BraKet[orbVec[i]].func_ptr->data.n1[0] != 0 and
                                BraKet[orbVec[j]].func_ptr->data.n1[0] != 0)
                                continue;
                            S(orbVec[i], orbVec[j]) += S_temp(i, j);
                        }
                    }
                }
            }
        }
        if (serial) {
#pragma omp critical
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < N; j++) { S(i, j) += S_omp(i, j); }
            }
        }
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j <= i; j++) {
            if (i != j) S(j, i) = std::conj(S(i, j)); // ensure exact symmetri
        }
    }

    // Assumes linearity: result is sum of all nodes contributions
    mrcpp::mpi::allreduce_matrix(S, mrcpp::mpi::comm_wrk);

    // multiply by CompFunction multiplicative factor
    ComplexVector Fac = ComplexVector::Zero(N);
    for (int i = 0; i < N; i++) {
        if (!mrcpp::mpi::my_func(BraKet[i])) continue;
        Fac[i] = BraKet[i].func_ptr->data.c1[0];
    }
    mrcpp::mpi::allreduce_vector(Fac, mrcpp::mpi::comm_wrk);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) { S(i, j) *= std::conj(Fac[i]) * Fac[j]; }
    }

    return S;
}

/** @brief Compute the overlap matrix S_ij = <bra_i|ket_j>
 *
 *  Will take the conjugate of bra before integrating
 */
ComplexMatrix calc_overlap_matrix_cplx(CompFunctionVector &Bra, CompFunctionVector &Ket) {
    mrcpp::mpi::barrier(mrcpp::mpi::comm_wrk); // for consistent timings
    bool braisreal = !Bra[0].iscomplex();
    bool ketisreal = !Ket[0].iscomplex();
    if (braisreal or ketisreal) {
        // temporary solution: copy as complex trees
        if (braisreal) {
            for (int i = 0; i < Bra.size(); i++) {
                Bra[i].CompD[0]->CopyTreeToComplex(Bra[i].CompC[0]);
                Bra[i].func_ptr->iscomplex = 1;
            }
        }
        if (ketisreal) {
            for (int i = 0; i < Ket.size(); i++) {
                Ket[i].CompD[0]->CopyTreeToComplex(Ket[i].CompC[0]);
                Ket[i].func_ptr->iscomplex = 1;
            }
        }
    }
    MultiResolutionAnalysis<3> *mra = Bra.vecMRA;

    int N = Bra.size();
    int M = Ket.size();
    ComplexMatrix S = ComplexMatrix::Zero(N, M);

    IntVector conjMatBra = IntVector::Zero(N);
    for (int i = 0; i < N; i++) {
        if (!mrcpp::mpi::my_func(Bra[i])) continue;
        conjMatBra[i] = (Bra[i].conjugate()) ? 1 : 0;
    }
    mrcpp::mpi::allreduce_vector(conjMatBra, mrcpp::mpi::comm_wrk);
    IntVector conjMatKet = IntVector::Zero(M);
    for (int i = 0; i < M; i++) {
        if (!mrcpp::mpi::my_func(Ket[i])) continue;
        conjMatKet[i] = (Ket[i].conjugate()) ? 1 : 0;
    }
    mrcpp::mpi::allreduce_vector(conjMatKet, mrcpp::mpi::comm_wrk);

    // 1) make union tree without coefficients for Bra (supposed smallest)
    mrcpp::FunctionTree<3> refTree(*mra);
    mrcpp::mpi::allreduce_Tree_noCoeff(refTree, Bra, mpi::comm_wrk);
    // note that Ket is not part of union grid: if a node is in ket but not in Bra, the dot product is zero.

    int sizecoeff = (1 << refTree.getDim()) * refTree.getKp1_d();
    int sizecoeffW = ((1 << refTree.getDim()) - 1) * refTree.getKp1_d();

    // get a list of all nodes in union grid, as defined by their indices
    std::vector<double *> coeffVec_ref;
    std::vector<int> indexVec_ref;    // serialIx of the nodes
    std::vector<int> parindexVec_ref; // serialIx of the parent nodes
    std::vector<double> scalefac;
    int max_ix;

    refTree.makeCoeffVector(coeffVec_ref, indexVec_ref, parindexVec_ref, scalefac, max_ix, refTree);
    int max_n = indexVec_ref.size();
    max_ix++;

    bool serial = mrcpp::mpi::wrk_size == 1; // flag for serial/MPI switch

    // only used for serial case:
    std::vector<std::vector<ComplexDouble *>> coeffVecBra(N);
    std::map<int, std::vector<int>> node2orbVecBra; // for each node index, gives a vector with the indices of the orbitals using this node
    std::vector<std::map<int, int>> orb2nodeBra(N); // for a given orbital and a given node, gives the node index in
                                                    // the orbital given the node index in the reference tree
    std::vector<std::vector<ComplexDouble *>> coeffVecKet(M);
    std::map<int, std::vector<int>> node2orbVecKet; // for each node index, gives a vector with the indices of the orbitals using this node
    std::vector<std::map<int, int>> orb2nodeKet(M); // for a given orbital and a given node, gives the node index in
                                                    // the orbital given the node index in the reference tree
    mrcpp::BankAccount nodesBra;
    mrcpp::BankAccount nodesKet;

    // In the serial case we store the coeff pointers in coeffVec. In the mpi case the coeff are stored in the bank
    if (serial) {
        // 2) make list of all coefficients, and their reference indices
        // for different orbitals, indexVec will give the same index for the same node in space
        // TODO? : do not copy coefficients, but use directly the pointers
        // could OMP parallelize, but is fast anyway
        std::vector<int> parindexVec; // serialIx of the parent nodes
        std::vector<int> indexVec;    // serialIx of the nodes
        for (int j = 0; j < N; j++) {
            // make vector with all coef pointers and their indices in the union grid
            Bra[j].complex().makeCoeffVector(coeffVecBra[j], indexVec, parindexVec, scalefac, max_ix, refTree);
            // make a map that gives j from indexVec
            int orb_node_ix = 0;
            for (int ix : indexVec) {
                orb2nodeBra[j][ix] = orb_node_ix++;
                if (ix < 0) continue;
                node2orbVecBra[ix].push_back(j);
            }
        }
        for (int j = 0; j < M; j++) {
            Ket[j].complex().makeCoeffVector(coeffVecKet[j], indexVec, parindexVec, scalefac, max_ix, refTree);
            // make a map that gives j from indexVec
            int orb_node_ix = 0;
            for (int ix : indexVec) {
                orb2nodeKet[j][ix] = orb_node_ix++;
                if (ix < 0) continue;
                node2orbVecKet[ix].push_back(j);
            }
        }

    } else { // MPI case
        // 2) send own nodes to bank, identifying them through the serialIx of refTree
        save_nodes(Bra, refTree, nodesBra);
        save_nodes(Ket, refTree, nodesKet);
        mrcpp::mpi::barrier(mrcpp::mpi::comm_wrk); // wait until everything is stored before fetching!
    }

    // 3) make dot product for all the nodes and accumulate into S
    int totsiz = 0;
    int totget = 0;
    int mxtotsiz = 0;
    int ibank = 0;
    // the omp crashes sometime for unknown reasons?
#pragma omp parallel if (serial)
    {
        ComplexMatrix S_omp = ComplexMatrix::Zero(N, M); // copy for each thread

#pragma omp for schedule(dynamic)
        for (int n = 0; n < max_n; n++) {
            if (n % mrcpp::mpi::wrk_size != mrcpp::mpi::wrk_rank) continue;
            int csize;
            std::vector<int> orbVecBra; // identifies which Bra orbitals use this node
            std::vector<int> orbVecKet; // identifies which Ket orbitals use this node
            if (parindexVec_ref[n] < 0)
                csize = sizecoeff;
            else
                csize = sizecoeffW;
            if (serial) {
                int node_ix = indexVec_ref[n];      // SerialIx for this node in the reference tree
                int shift = sizecoeff - sizecoeffW; // to copy only wavelet part
                ComplexMatrix coeffBlockBra(csize, node2orbVecBra[node_ix].size());
                ComplexMatrix coeffBlockKet(csize, node2orbVecKet[node_ix].size());
                if (parindexVec_ref[n] < 0) shift = 0;

                for (int j : node2orbVecBra[node_ix]) { // loop over indices of the orbitals using this node
                    int orb_node_ix = orb2nodeBra[j][node_ix];
                    for (int k = 0; k < csize; k++) coeffBlockBra(k, orbVecBra.size()) = coeffVecBra[j][orb_node_ix][k + shift];
                    orbVecBra.push_back(j);
                }
                for (int j : node2orbVecKet[node_ix]) { // loop over indices of the orbitals using this node
                    int orb_node_ix = orb2nodeKet[j][node_ix];
                    for (int k = 0; k < csize; k++) coeffBlockKet(k, orbVecKet.size()) = coeffVecKet[j][orb_node_ix][k + shift];
                    orbVecKet.push_back(j);
                }

                if (orbVecBra.size() > 0 and orbVecKet.size() > 0) {
                    ComplexMatrix S_temp(orbVecBra.size(), orbVecKet.size());
                    if (not conjMatBra[0] and not conjMatBra[0]) {
                        S_temp.noalias() = coeffBlockBra.transpose().conjugate() * coeffBlockKet;
                    } else if (conjMatBra[0] and not conjMatBra[0]) {
                        S_temp.noalias() = coeffBlockBra.transpose() * coeffBlockKet;
                    } else if (not conjMatBra[0] and conjMatBra[0]) {
                        S_temp.noalias() = coeffBlockBra.transpose() * coeffBlockKet.transpose();
                    } else if (conjMatBra[0] and conjMatBra[0]) {
                        S_temp.noalias() = coeffBlockBra * coeffBlockKet.transpose();
                    } else
                        MSG_ABORT("Unexpected case");
                    for (int i = 0; i < orbVecBra.size(); i++) {
                        for (int j = 0; j < orbVecKet.size(); j++) {
                            if (Bra[orbVecBra[i]].func_ptr->data.n1[0] != Ket[orbVecKet[j]].func_ptr->data.n1[0] and Bra[orbVecBra[i]].func_ptr->data.n1[0] != 0 and
                                Ket[orbVecKet[j]].func_ptr->data.n1[0] != 0)
                                continue;
                            S_omp(orbVecBra[i], orbVecKet[j]) += S_temp(i, j);
                        }
                    }
                }
            } else { // MPI case

                ComplexMatrix coeffBlockBra(csize, N);
                ComplexMatrix coeffBlockKet(csize, M);
                nodesBra.get_nodeblock(indexVec_ref[n], coeffBlockBra.data(), orbVecBra); // get Bra parts
                nodesKet.get_nodeblock(indexVec_ref[n], coeffBlockKet.data(), orbVecKet); // get Ket parts
                totsiz += orbVecBra.size() * orbVecKet.size();
                mxtotsiz += N * M;
                totget += orbVecBra.size() + orbVecKet.size();
                if (orbVecBra.size() > 0 and orbVecKet.size() > 0) {
                    ComplexMatrix S_temp(orbVecBra.size(), orbVecKet.size());
                    coeffBlockBra.conservativeResize(Eigen::NoChange, orbVecBra.size());
                    coeffBlockKet.conservativeResize(Eigen::NoChange, orbVecKet.size());
                    if (not conjMatBra[0] and not conjMatBra[0]) {
                        S_temp.noalias() = coeffBlockBra.transpose().conjugate() * coeffBlockKet;
                    } else if (conjMatBra[0] and not conjMatBra[0]) {
                        S_temp.noalias() = coeffBlockBra.transpose() * coeffBlockKet;
                    } else if (not conjMatBra[0] and conjMatBra[0]) {
                        S_temp.noalias() = coeffBlockBra.transpose() * coeffBlockKet.transpose();
                    } else if (conjMatBra[0] and conjMatBra[0]) {
                        S_temp.noalias() = coeffBlockBra * coeffBlockKet.transpose();
                    } else
                        MSG_ABORT("Unexpected case");

                    for (int i = 0; i < orbVecBra.size(); i++) {
                        for (int j = 0; j < orbVecKet.size(); j++) {
                            if (Bra[orbVecBra[i]].func_ptr->data.n1[0] != Ket[orbVecKet[j]].func_ptr->data.n1[0] and Bra[orbVecBra[i]].func_ptr->data.n1[0] != 0 and
                                Ket[orbVecKet[j]].func_ptr->data.n1[0] != 0)
                                continue;
                            S(orbVecBra[i], orbVecKet[j]) += S_temp(i, j);
                        }
                    }
                }
            }
        }
        if (serial) {
#pragma omp critical
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < M; j++) { S(i, j) += S_omp(i, j); }
            }
        }
    }

    // 4) collect results from all MPI. Linearity: result is sum of all node contributions

    mrcpp::mpi::allreduce_matrix(S, mrcpp::mpi::comm_wrk);

    // multiply by CompFunction multiplicative factor
    ComplexVector FacBra = ComplexVector::Zero(N);
    ComplexVector FacKet = ComplexVector::Zero(M);
    for (int i = 0; i < N; i++) {
        if (!mrcpp::mpi::my_func(Bra[i])) continue;
        FacBra[i] = Bra[i].func_ptr->data.c1[0];
    }
    for (int i = 0; i < M; i++) {
        if (!mrcpp::mpi::my_func(Ket[i])) continue;
        FacKet[i] = Ket[i].func_ptr->data.c1[0];
    }
    mrcpp::mpi::allreduce_vector(FacBra, mrcpp::mpi::comm_wrk);
    mrcpp::mpi::allreduce_vector(FacKet, mrcpp::mpi::comm_wrk);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) { S(i, j) *= std::conj(FacBra[i]) * FacKet[j]; }
    }

    // restore input
    if (braisreal) {
        for (int i = 0; i < Bra.size(); i++) {
            delete Bra[i].CompC[0];
            Bra[i].CompC[0] = nullptr;
            Bra[i].func_ptr->iscomplex = 0;
            Bra[i].func_ptr->isreal = 1;
        }
    }
    if (ketisreal) {
        for (int i = 0; i < Ket.size(); i++) {
            delete Ket[i].CompC[0];
            Ket[i].CompC[0] = nullptr;
            Ket[i].func_ptr->iscomplex = 0;
            Ket[i].func_ptr->isreal = 1;
        }
    }
    return S;
}

/** @brief Compute the overlap matrix S_ij = <bra_i|ket_j>
 *
 */
ComplexMatrix calc_overlap_matrix(CompFunctionVector &Bra, CompFunctionVector &Ket) {

    if (Bra[0].iscomplex() or Ket[0].iscomplex()) { return calc_overlap_matrix_cplx(Bra, Ket); }

    mrcpp::mpi::barrier(mrcpp::mpi::comm_wrk); // for consistent timings

    MultiResolutionAnalysis<3> *mra = Bra.vecMRA;

    int N = Bra.size();
    int M = Ket.size();
    ComplexMatrix S = ComplexMatrix::Zero(N, M);

    // 1) make union tree without coefficients for Bra (supposed smallest)
    mrcpp::FunctionTree<3> refTree(*mra);
    mrcpp::mpi::allreduce_Tree_noCoeff(refTree, Bra, mpi::comm_wrk);
    // note that Ket is not part of union grid: if a node is in ket but not in Bra, the dot product is zero.

    int sizecoeff = (1 << refTree.getDim()) * refTree.getKp1_d();
    int sizecoeffW = ((1 << refTree.getDim()) - 1) * refTree.getKp1_d();

    // get a list of all nodes in union grid, as defined by their indices
    std::vector<double *> coeffVec_ref;
    std::vector<int> indexVec_ref;    // serialIx of the nodes
    std::vector<int> parindexVec_ref; // serialIx of the parent nodes
    std::vector<double> scalefac;
    int max_ix;

    refTree.makeCoeffVector(coeffVec_ref, indexVec_ref, parindexVec_ref, scalefac, max_ix, refTree);
    int max_n = indexVec_ref.size();
    max_ix++;

    bool serial = mrcpp::mpi::wrk_size == 1; // flag for serial/MPI switch

    // only used for serial case:
    std::vector<std::vector<double *>> coeffVecBra(N);
    std::map<int, std::vector<int>> node2orbVecBra; // for each node index, gives a vector with the indices of the orbitals using this node
    std::vector<std::map<int, int>> orb2nodeBra(N); // for a given orbital and a given node, gives the node index in
                                                    // the orbital given the node index in the reference tree
    std::vector<std::vector<double *>> coeffVecKet(M);
    std::map<int, std::vector<int>> node2orbVecKet; // for each node index, gives a vector with the indices of the orbitals using this node
    std::vector<std::map<int, int>> orb2nodeKet(M); // for a given orbital and a given node, gives the node index in
                                                    // the orbital given the node index in the reference tree
    mrcpp::BankAccount nodesBra;
    mrcpp::BankAccount nodesKet;
    // In the serial case we store the coeff pointers in coeffVec. In the mpi case the coeff are stored in the bank
    if (serial) {
        // 2) make list of all coefficients, and their reference indices
        // for different orbitals, indexVec will give the same index for the same node in space
        // TODO? : do not copy coefficients, but use directly the pointers
        // could OMP parallelize, but is fast anyway
        std::vector<int> parindexVec; // serialIx of the parent nodes
        std::vector<int> indexVec;    // serialIx of the nodes
        for (int j = 0; j < N; j++) {
            // make vector with all coef pointers and their indices in the union grid
            Bra[j].real().makeCoeffVector(coeffVecBra[j], indexVec, parindexVec, scalefac, max_ix, refTree);
            // make a map that gives j from indexVec
            int orb_node_ix = 0;
            for (int ix : indexVec) {
                orb2nodeBra[j][ix] = orb_node_ix++;
                if (ix < 0) continue;
                node2orbVecBra[ix].push_back(j);
            }
        }
        for (int j = 0; j < M; j++) {
            Ket[j].real().makeCoeffVector(coeffVecKet[j], indexVec, parindexVec, scalefac, max_ix, refTree);
            // make a map that gives j from indexVec
            int orb_node_ix = 0;
            for (int ix : indexVec) {
                orb2nodeKet[j][ix] = orb_node_ix++;
                if (ix < 0) continue;
                node2orbVecKet[ix].push_back(j);
            }
        }

    } else { // MPI case
        // 2) send own nodes to bank, identifying them through the serialIx of refTree
        save_nodes(Bra, refTree, nodesBra);
        save_nodes(Ket, refTree, nodesKet);
        mrcpp::mpi::barrier(mrcpp::mpi::comm_wrk); // wait until everything is stored before fetching!
    }

    // 3) make dot product for all the nodes and accumulate into S
    int totsiz = 0;
    int totget = 0;
    int mxtotsiz = 0;
    int ibank = 0;
#pragma omp parallel if (serial)
    {
        DoubleMatrix S_omp = DoubleMatrix::Zero(N, M); // copy for each thread
        // NB: dynamic does give strange errors?
#pragma omp for schedule(static)
        for (int n = 0; n < max_n; n++) {
            if (n % mrcpp::mpi::wrk_size != mrcpp::mpi::wrk_rank) continue;
            int csize;
            std::vector<int> orbVecBra; // identifies which Bra orbitals use this node
            std::vector<int> orbVecKet; // identifies which Ket orbitals use this node
            if (parindexVec_ref[n] < 0)
                csize = sizecoeff;
            else
                csize = sizecoeffW;
            if (serial) {
                int node_ix = indexVec_ref[n];      // SerialIx for this node in the reference tree
                int shift = sizecoeff - sizecoeffW; // to copy only wavelet part
                DoubleMatrix coeffBlockBra(csize, node2orbVecBra[node_ix].size());
                DoubleMatrix coeffBlockKet(csize, node2orbVecKet[node_ix].size());
                if (parindexVec_ref[n] < 0) shift = 0;

                for (int j : node2orbVecBra[node_ix]) { // loop over indices of the orbitals using this node
                    int orb_node_ix = orb2nodeBra[j][node_ix];
                    for (int k = 0; k < csize; k++) coeffBlockBra(k, orbVecBra.size()) = coeffVecBra[j][orb_node_ix][k + shift];
                    orbVecBra.push_back(j);
                }
                for (int j : node2orbVecKet[node_ix]) { // loop over indices of the orbitals using this node
                    int orb_node_ix = orb2nodeKet[j][node_ix];
                    for (int k = 0; k < csize; k++) coeffBlockKet(k, orbVecKet.size()) = coeffVecKet[j][orb_node_ix][k + shift];
                    orbVecKet.push_back(j);
                }

                if (orbVecBra.size() > 0 and orbVecKet.size() > 0) {
                    DoubleMatrix S_temp(orbVecBra.size(), orbVecKet.size());
                    S_temp.noalias() = coeffBlockBra.transpose() * coeffBlockKet;

                    for (int i = 0; i < orbVecBra.size(); i++) {
                        for (int j = 0; j < orbVecKet.size(); j++) {
                            if (Bra[orbVecBra[i]].func_ptr->data.n1[0] != Ket[orbVecKet[j]].func_ptr->data.n1[0] and Bra[orbVecBra[i]].func_ptr->data.n1[0] != 0 and
                                Ket[orbVecKet[j]].func_ptr->data.n1[0] != 0)
                                continue;
                            S_omp(orbVecBra[i], orbVecKet[j]) += S_temp(i, j);
                        }
                    }
                }
            } else { // MPI case

                DoubleMatrix coeffBlockBra(csize, N);
                DoubleMatrix coeffBlockKet(csize, M);
                nodesBra.get_nodeblock(indexVec_ref[n], coeffBlockBra.data(), orbVecBra); // get Bra parts
                nodesKet.get_nodeblock(indexVec_ref[n], coeffBlockKet.data(), orbVecKet); // get Ket parts
                totsiz += orbVecBra.size() * orbVecKet.size();
                mxtotsiz += N * M;
                totget += orbVecBra.size() + orbVecKet.size();
                if (orbVecBra.size() > 0 and orbVecKet.size() > 0) {
                    DoubleMatrix S_temp(orbVecBra.size(), orbVecKet.size());
                    coeffBlockBra.conservativeResize(Eigen::NoChange, orbVecBra.size());
                    coeffBlockKet.conservativeResize(Eigen::NoChange, orbVecKet.size());
                    S_temp.noalias() = coeffBlockBra.transpose() * coeffBlockKet;
                    for (int i = 0; i < orbVecBra.size(); i++) {
                        for (int j = 0; j < orbVecKet.size(); j++) {
                            if (Bra[orbVecBra[i]].func_ptr->data.n1[0] != Ket[orbVecKet[j]].func_ptr->data.n1[0] and Bra[orbVecBra[i]].func_ptr->data.n1[0] != 0 and
                                Ket[orbVecKet[j]].func_ptr->data.n1[0] != 0)
                                continue;
                            S(orbVecBra[i], orbVecKet[j]) += S_temp(i, j);
                        }
                    }
                }
            }
        }
        if (serial) {
#pragma omp critical
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < M; j++) { S(i, j) += S_omp(i, j); }
            }
        }
    }

    // 4) collect results from all MPI. Linearity: result is sum of all node contributions

    mrcpp::mpi::allreduce_matrix(S, mrcpp::mpi::comm_wrk);

    // multiply by CompFunction multiplicative factor
    ComplexVector FacBra = ComplexVector::Zero(N);
    ComplexVector FacKet = ComplexVector::Zero(M);
    for (int i = 0; i < N; i++) {
        if (!mrcpp::mpi::my_func(Bra[i])) continue;
        FacBra[i] = Bra[i].func_ptr->data.c1[0];
    }
    for (int i = 0; i < M; i++) {
        if (!mrcpp::mpi::my_func(Ket[i])) continue;
        FacKet[i] = Ket[i].func_ptr->data.c1[0];
    }
    mrcpp::mpi::allreduce_vector(FacBra, mrcpp::mpi::comm_wrk);
    mrcpp::mpi::allreduce_vector(FacKet, mrcpp::mpi::comm_wrk);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) { S(i, j) *= std::conj(FacBra[i]) * FacKet[j]; }
    }

    return S;
}

/** @brief Orthogonalize the functions in Bra against all orbitals in Ket
 *
 */
void orthogonalize(double prec, CompFunctionVector &Bra, CompFunctionVector &Ket) {
    // TODO: generalize for cases where Ket functions are not orthogonal to each other?
    ComplexMatrix S = calc_overlap_matrix(Bra, Ket);
    int N = Bra.size();
    int M = Ket.size();
    DoubleVector Ketnorms = DoubleVector::Zero(M);
    for (int i = 0; i < M; i++) {
        if (mpi::my_func(Ket[i])) Ketnorms(i) = Ket[i].getSquareNorm();
    }
    mrcpp::mpi::allreduce_vector(Ketnorms, mrcpp::mpi::comm_wrk);
    ComplexMatrix rmat = ComplexMatrix::Zero(M, N);
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < M; i++) { rmat(i, j) = 0.0 - S.conjugate()(j, i) / Ketnorms(i); }
    }
    CompFunctionVector rotatedKet(N);
    rotate(Ket, rmat, rotatedKet, prec / M);
    for (int j = 0; j < N; j++) {
        if (my_func(Bra[j])) Bra[j].add(1.0, rotatedKet[j]);
    }
}

/** @brief Orthogonalize the Bra against Ket
 *
 */
template <int D> void orthogonalize(double prec, CompFunction<D> &Bra, CompFunction<D> &Ket) {
    ComplexDouble overlap = dot(Bra, Ket);
    double sq_norm = Ket.getSquareNorm();
    for (int i = 0; i < Bra.Ncomp(); i++) {
        if (Bra.isreal()) {
            if (abs(overlap.imag()) > MachineZero) MSG_ABORT("NOT IMPLEMENTED");
            Bra.CompD[i]->add_inplace(-overlap.real() / sq_norm, *Ket.CompD[i]);
        } else {
            if (Ket.isreal()) MSG_ABORT("NOT IMPLEMENTED");
            Bra.CompC[i]->add_inplace(-std::conj(overlap / sq_norm), *Ket.CompC[i]);
            overlap = dot(Bra, Ket);
        }
    }
}

template ComplexDouble dot(CompFunction<3> bra, CompFunction<3> ket);
template void project(CompFunction<3> &out, RepresentableFunction<3, double> &f, double prec);
template void project(CompFunction<3> &out, RepresentableFunction<3, ComplexDouble> &f, double prec);
template void multiply(CompFunction<3> &out, CompFunction<3> inp_a, CompFunction<3> inp_b, double prec, bool absPrec, bool useMaxNorms, bool conjugate);
template void multiply(CompFunction<3> &out, FunctionTree<3, double> &inp_a, RepresentableFunction<3, double> &f, double prec, int nrefine = 0, bool conjugate);
template void multiply(CompFunction<3> &out, FunctionTree<3, ComplexDouble> &inp_a, RepresentableFunction<3, ComplexDouble> &f, double prec, int nrefine = 0, bool conjugate);
template void multiply(CompFunction<3> &out, CompFunction<3> &inp_a, RepresentableFunction<3, double> &f, double prec, int nrefine = 0, bool conjugate);
template void multiply(CompFunction<3> &out, CompFunction<3> &inp_a, RepresentableFunction<3, ComplexDouble> &f, double prec, int nrefine = 0, bool conjugate);
template void deep_copy(CompFunction<3> *out, const CompFunction<3> &inp);
template void deep_copy(CompFunction<3> &out, const CompFunction<3> &inp);
template void add(CompFunction<3> &out, ComplexDouble a, CompFunction<3> inp_a, ComplexDouble b, CompFunction<3> inp_b, double prec, bool conjugate);
template void linear_combination(CompFunction<3> &out, const std::vector<ComplexDouble> &c, std::vector<CompFunction<3>> &inp, double prec, bool conjugate);
template double node_norm_dot(CompFunction<3> bra, CompFunction<3> ket);
template void orthogonalize(double prec, CompFunction<3> &Bra, CompFunction<3> &Ket);

} // namespace mrcpp
