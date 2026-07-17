/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2021 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
 *
 * This file is part of MRCPP.
 *
 * MRCPP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRCPP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRCPP.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRCPP, see:
 * <https://mrcpp.readthedocs.io/>
 */

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

#include <any>

/* Some rules for CompFunction:
 * NComp is the number of components. If Ncomp>0, the corresponding trees must exist (can be only empty roots).
 * The other trees should be set to nullptr.
 * The trees and data can be shared among several CompFunction; this is managed automatically by "std::make_shared"
 * Normally the CompFunction must be eiher real or complex (or none if none is defined anyway).
 * Though it is allowed in some cases to have both and the code should preferably allow this. (It is used temporary
 * when we need a Complex type, but the trees are real: the tree is then copied as a complex tree in the same CompFunction).
 * TreePtr (aka func_ptr) is the part potentially shared with others with "std::make_shared". It contains the pointers to the trees.
 * The static data (number of components, real/complex, conjugaison, integers used for spin etc.) are store in func_ptr.data.
 */

namespace mrcpp {

// template <int D> MultiResolutionAnalysis<D> *defaultCompMRA = nullptr; // Global MRA RAW_VER
template <int D> std::shared_ptr<MultiResolutionAnalysis<D>> defaultCompMRA; // Global MRA would prolly be better but requires many changes

template <int D> CompFunction<D>::CompFunction(MultiResolutionAnalysis<D> &mra) {
    defaultCompMRA<D> = std::make_shared<MultiResolutionAnalysis<D>>(mra);
    func_ptr = std::make_shared<TreePtr<D>>(false);
    CompD = func_ptr->real;
    CompC = func_ptr->cplx;
    for (int i = 0; i < 4; i++) CompD[i] = nullptr;
    for (int i = 0; i < 4; i++) CompC[i] = nullptr;
}

template <int D> CompFunction<D>::CompFunction(MultiResolutionAnalysis<D> &mra, int nComponents) {
    defaultCompMRA<D> = std::make_shared<MultiResolutionAnalysis<D>>(mra);
    func_ptr = std::make_shared<TreePtr<D>>(false);
    CompD = func_ptr->real;
    CompC = func_ptr->cplx;
    func_ptr->data.Ncomp = nComponents;
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
 * @param n1: for instance spin (1=alpha, -1=beta, 2=paired) or anything else
 */
template <int D> CompFunction<D>::CompFunction(int n1, int nComponents) {
    // std::cout << "CompFunction::CompFunction constructor with n1: start " << n1 << std::endl;
    func_ptr = std::make_shared<TreePtr<D>>(false);
    CompD = func_ptr->real;
    CompC = func_ptr->cplx;
    for (int i = 0; i < 4; i++) CompD[i] = nullptr;
    for (int i = 0; i < 4; i++) CompC[i] = nullptr;
    func_ptr->data.n1[0] = n1; // for orbitals: spin 
    func_ptr->data.n2[0] = -1; 
    func_ptr->data.n3[0] = 0;
    func_ptr->rank = 0;
    func_ptr->isreal = 1;
    func_ptr->iscomplex = 0;
    func_ptr->data.shared = false;
    // std::cout << "CompFunction::CompFunction constructor with n1: " << n1 << std::endl;
    func_ptr->data.Ncomp = nComponents;
}

/*
 * Empty functions (no components defined)
 * Not universal like the other constructors, but specific to set up spinors
 * @param spin: "Paired", "Alpha", "Beta"
 * @param nComponents: number of components (1 for scalar, 2 for spinor, 4 for Dirac spinor)
 * @param spin_type: "Large", "Small"
 */
// template <int D> CompFunction<D>::CompFunction(std::string spin, int nComponents, std::string spin_type) {
//     // --- Setting up the pointers to the trees
//     func_ptr = std::make_shared<TreePtr<D>>(false);
//     CompD = func_ptr->real;
//     CompC = func_ptr->cplx;
//     for (int i = 0; i < 4; i++) CompD[i] = nullptr;
//     for (int i = 0; i < 4; i++) CompC[i] = nullptr;

//     // --- Setting up spin
//     // first determine if large or small component
//     int spin_index = 0;
//     if (spin_type == "Large" || spin_type == "large" || spin_type == "L" || spin_type == "l") { // Default, Keeping this here (even though it's redundant) just to make the code more readable
//         spin_index = 0;
//     }
//     else if (spin_type == "Small" || spin_type == "small" || spin_type == "S" || spin_type == "s") {
//         spin_index = 2;
//     }
//     else {
//         MSG_ERROR( "CompFunction: unknown component type ");
//     }
//     if (spin == "Paired" || spin == "paired" || spin == "P" || spin == "p")
//         func_ptr->data.n1[spin_index] = 2;
//     else if (spin == "Alpha" || spin == "alpha" || spin == "A" || spin == "a")
//         func_ptr->data.n1[spin_index] = 1;
//     else if (spin == "Beta" || spin == "beta" || spin == "B" || spin == "b")
//         func_ptr->data.n1[spin_index + 1] = 1;
//     else {
//         MSG_ERROR( "CompFunction: unknown spin type ");
//     }
//     // func_ptr->data.n1[0] = n1;
//     // func_ptr->data.n2[0] = -1;
//     // func_ptr->data.n3[0] = 0;
//     func_ptr->rank = 0;
//     func_ptr->isreal = 1;
//     func_ptr->iscomplex = 0;
//     func_ptr->data.shared = false;
//     func_ptr->data.Ncomp = nComponents;
// }

/*
 * Empty functions (no components defined)
 */
template <int D> CompFunction<D>::CompFunction(int n1, bool share, int nComponents) {
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
    func_ptr->data.Ncomp = nComponents;
    // for (int i = 0; i < 4; i++) CompD[i]->setMRA() = *defaultCompMRA<D>;
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
    func_ptr->data.Ncomp = compfunc.func_ptr->data.Ncomp;
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
        func_ptr->data.Ncomp = compfunc.func_ptr->data.Ncomp;
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

template <int D> ComplexDouble CompFunction<D>::integrateSide(int dim, bool positiveSide) const {
    ComplexDouble integral;
    if (isreal())
        integral = CompD[0]->integrateSide(dim, positiveSide);
    else
        integral = CompC[0]->integrateSide(dim, positiveSide);
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

template <int D> void CompFunction<D>::calcSquareNorm() {
    for (int i = 0; i < Ncomp(); i++) {
        if (isreal() and CompD[i] != nullptr) {
            CompD[i]->calcSquareNorm();
        } else if (iscomplex() and CompC[i] != nullptr) {
            CompC[i]->calcSquareNorm();
        }
    }
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

//  @brief Allocate one empty trees for one specific component.
//  The tree must be defined as real or complex already.
//  @param ialloc is index allocated. ialloc=0 allocates the tree with index zero.
//  deletes the old tree if found.
template <int D> void CompFunction<D>::alloc_comp(int ialloc, bool zero) {
    if (defaultCompMRA<D> == nullptr) MSG_ABORT("Default MRA not yet defined");
    if (isreal() == 0 and iscomplex() == 0) MSG_ABORT("Function must be defined either real or complex");
    int i = ialloc;
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

// template <int D> void CompFunction<D>::freeTrees()

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
    // return this->complex()
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
    CompD[i] = tree;
    if (tree != nullptr) {
        func_ptr->Ncomp = std::max(Ncomp(), i + 1);
    } else {
        func_ptr->Ncomp = std::min(Ncomp(), i);
    }
}

template <int D> void CompFunction<D>::setCplx(FunctionTree<D, ComplexDouble> *tree, int i) {
    func_ptr->iscomplex = 1;
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
    //check if the prefactors are zero, in which case,
    //delete tree, reallocate to zero and reset prefactor to one
    for (int i = 0; i < inp.Ncomp(); i++) {
        if (std::abs((*this).func_ptr->data.c1[i].real())<MachineZero and std::abs((*this).func_ptr->data.c1[i].imag()) < MachineZero) {
            this->alloc_comp(i, true);
            (*this).func_ptr->data.c1[i] = {1.0,0.0};
        }
    }
    if (inp.getSquareNorm() < MachineZero) {
        // nothing to add
    } else if((this->getSquareNorm() < MachineZero) or (Ncomp() < inp.Ncomp())) {
        //self empty, copy inp into self
        func_ptr->data = inp.func_ptr->data;
        alloc(inp.Ncomp(), true);
        for (int i = 0; i < inp.Ncomp(); i++) {
            //equivalent to copying and rescaling
            //No need to include the prefactors c1 in this addition
            //because they have been inherited from the inp
            if(inp.isreal() and (std::abs(c.imag())<MachineZero)) CompD[i]->add_inplace(c.real(), *inp.CompD[i]);
            else {
                if (inp.isreal()) { //provision in case inp is real and the coefficient c is not
                    inp.CompC[i] = inp.CompD[i]->CopyTreeToComplex();
                    //out's type is copied from inp
                    (*this).CompC[i] = (*this).CompD[i]->CopyTreeToComplex();
                    delete (*this).CompD[i];
                    (*this).CompD[i] = nullptr;
                    this->defcomplex(); //will be called up to 4 times, but that's ok
                }
                CompC[i]->add_inplace(c, *inp.CompC[i]);
                if (inp.isreal()){ //restoring inp to what is was
                    delete inp.CompC[i];
                    inp.CompC[i] = nullptr;
                }
            } 
        }
        //Set the norm to be what it should be, instead of the default 0
        this->calcSquareNorm();
    } else {
        //adding inp to non-empty self
        for (int i = 0; i < inp.Ncomp(); i++) {
            if (this->isreal() and inp.isreal() and std::abs((c*inp.func_ptr->data.c1[i]/(*this).func_ptr->data.c1[i]).imag()) < MachineZero) {
                //inplace addition, thus the output keeps its prefactor c1
                //c1*out <- c1*out + c2*inp  means c1*out = c1*(out + c2*inp/c1)
                CompD[i]->add_inplace((c*inp.func_ptr->data.c1[i]/(*this).func_ptr->data.c1[i]).real(), *inp.CompD[i]);
            } else {
                if (this->isreal()) {
                    for (int comp = 0; comp < this->Ncomp(); comp++){
                        CompC[comp] = CompD[comp]->CopyTreeToComplex();
                        delete CompD[comp];
                        CompD[comp] = nullptr;
                    }
                    func_ptr->iscomplex = true;
                    func_ptr->isreal = false;
                }
                if (inp.isreal()) { //provision in case inp is real and the coefficient c is not
                    inp.CompC[i] = inp.CompD[i]->CopyTreeToComplex();
                }
                //as explained in the real case:
                //c1*out <- c1*out + c2*inp  means c1*out = c1*(out + c2*inp/c1)
                CompC[i]->add_inplace(c*inp.func_ptr->data.c1[i]/(*this).func_ptr->data.c1[i], *inp.CompC[i]);
                if (inp.isreal()){ //restoring inp to what is was
                    delete inp.CompC[i];
                    inp.CompC[i] = nullptr;
                }
            }
        }
    }
}

template <int D> void CompFunction<D>::upgradeToComplex(){
    for (int i = 0; i < Ncomp(); i++) {
        CompC[i] = CompD[i]->CopyTreeToComplex();
        delete CompD[i];
        CompD[i] = nullptr;
    }
    func_ptr->iscomplex = 1;
    func_ptr->isreal = 0;
}    
template <int D> int CompFunction<D>::crop(double prec, bool absPrec) {
    if (prec < 0.0) return 0;
    int nChunksremoved = 0;
    for (int i = 0; i < Ncomp(); i++) {
        if (isreal()) {
            nChunksremoved += CompD[i]->crop(prec, 1.0, absPrec);
        } else {
            nChunksremoved += CompC[i]->crop(prec, 1.0, absPrec);
        }
    }
    return nChunksremoved;
}

/** @brief In place multiply with scalar. Fully in-place.*/
template <int D> void CompFunction<D>::rescale(ComplexDouble c) {
    bool need_to_rescale = not(isShared()) or mpi::share_master();
    bool need_complex_functions = ((abs(c.imag()) > MachineZero) and isreal());
    if(need_complex_functions) upgradeToComplex();
    if (need_to_rescale) {
        if (iscomplex()) {
            for (int i = 0; i < Ncomp(); i++) {
                CompC[i]->rescale(c);
            }
        } else {
            for (int i = 0; i < Ncomp(); i++) {
                if (abs(c.imag()) > MachineZero) { 
                    CompC[i] = CompD[i]->CopyTreeToComplex();
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
        // for (int i = 0; i < Ncomp(); i++) {
        //     if (iscomplex()) {
        //         CompC[i]->rescale(c);
        //     } else {
        //         if (abs(c.imag()) > MachineZero) { // works only only for NComp==1)
        //             CompD[i]->CopyTreeToComplex(CompC[i]);
        //             delete CompD[i];
        //             CompD[i] = nullptr;
        //             func_ptr->iscomplex = true;
        //             func_ptr->isreal = false;
        //             CompC[i]->rescale(c);
        //         } else {
        //             CompD[i]->rescale(c.real());
        //         }
        //     }
        // }
    } else
        MSG_ABORT("Not implemented");
}

template class MultiResolutionAnalysis<1>;
template class MultiResolutionAnalysis<2>;
template class MultiResolutionAnalysis<3>;
template class CompFunction<1>;
template class CompFunction<2>;
template class CompFunction<3>;

/** @brief Deep copy that changes type from real to complex
 *
 * Deep copy: makes an exact copy with type complex from a real input
 */
template <int D> void CopyToComplex(CompFunction<D> &out, const CompFunction<D> &inp) {
    out.func_ptr->data = inp.func_ptr->data;
    out.defcomplex();
    // out.func_ptr->data.isreal = 0;
    // out.alloc(inp.Ncomp()); //unnecessary and might cause a small memory leak
    if (inp.getNNodes() == 0) return;
    for (int i = 0; i < inp.Ncomp(); i++) {
        if (inp.isreal()) {
            out.CompC[i] = inp.CompD[i]->CopyTreeToComplex();
        } else {
            // out.CompC[i] = inp.CompC[i]->deep_copy();
            //we should uniformise the syntax of functions in this code
            inp.CompC[i]->deep_copy(out.CompC[i]); 
        }
    }
}

/** @brief Deep copy (pointer version)
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

/** @brief Deep copy (reference version)
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
    // The coefficients c1 of inp_a,b are not relevant in this function and will be handled in linear_combination
    std::vector<ComplexDouble> coefs(2);
    coefs[0] = a;
    coefs[1] = b;


    std::vector<CompFunction<D>> funcs; // NB: not a CompFunctionVector, because not run in parallel!
    // If one of the inputs also serves as the output, we need to deep copy it to avoid 
    // uncontrolled behaviour in the linear_combination function
    if (&out == &inp_a) { //useless as we have a shallow copy in this function, instead of the code standard passing by reference
        CompFunction<D> out_a;
        deep_copy(out_a, inp_a);
        funcs.push_back(out_a);
    } else {
        CompFunction<D> out_a; //debug test
        deep_copy(out_a, inp_a); //debug test
        funcs.push_back(out_a); //debug test
        // funcs.push_back(inp_a);
    }
    if (&out == &inp_b) {
        CompFunction<D> out_b;
        deep_copy(out_b, inp_b);
        funcs.push_back(out_b);
    } else {
        CompFunction<D> out_b; //debug test
        deep_copy(out_b, inp_b); //debug test
        funcs.push_back(out_b); // debug test
        // funcs.push_back(inp_b);
    }

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
        for (int comp = 0; comp < inp[i].Ncomp(); comp++)
            if (inp[i].iscomplex() or std::abs((c[i]*inp[i].func_ptr->data.c1[comp]).imag()) > MachineZero) iscomplex = true;
    if (iscomplex) {
        out.func_ptr->data.iscomplex = 1;
        out.func_ptr->data.isreal = 0;
    }
    out.alloc(inp[0].Ncomp());
    for (int comp = 0; comp < inp[0].Ncomp(); comp++) {
        if (not iscomplex) {
            FunctionTreeVector<D, double> fvec; // one component vector
            for (int i = 0; i < inp.size(); i++) {
                if (std::norm(c[i]*inp[i].func_ptr->data.c1[comp]) < thrs) continue;
                if (inp[i].getNNodes() == 0 or inp[i].CompD[comp]->getSquareNorm() < thrs) continue;
                fvec.push_back(std::make_tuple((c[i]*inp[i].func_ptr->data.c1[comp]).real(), inp[i].CompD[comp]));
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
                    //We need to change all components to complex if we define the CompFunction as complex.
                    for (int kcomp = 0; kcomp < inp[i].Ncomp(); kcomp++){
                        inp[i].CompC[kcomp] = inp[i].CompD[kcomp]->CopyTreeToComplex();
                        delete inp[i].CompD[kcomp];
                        inp[i].CompD[kcomp] = nullptr; 
                    }
                    inp[i].defcomplex();
                    inp[i].func_ptr->isreal = 0;
                }
                if (std::norm(c[i]*inp[i].func_ptr->data.c1[comp]) < thrs) continue;
                if (inp[i].getNNodes() == 0 or inp[i].CompC[comp]->getSquareNorm() < thrs) continue;
                fvec.push_back(std::make_tuple(c[i]*inp[i].func_ptr->data.c1[comp], inp[i].CompC[comp]));
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
        //resetting the prefactor as it has been inherited from inp[0] (and thus might be different than 1)
        //and has been accounted for in the linear combination coefficients
        out.func_ptr->data.c1[comp] = {1.0,0.0};
    }
}

/** @brief out = conj(inp) * inp
 * 
 *  @param contrib: 4-vector of boolean selecting which components contribute to the density. By default a 4-vector of true.
 * 
 *  Note that output is always real
 * 
 *
 */
template <int D> void make_density(CompFunction<D> &out, CompFunction<D> &inp, double prec, std::vector<bool> contrib) {
    //provision in case inadequate contribution vector is provided, we just pad it with false
    if (contrib.size() < inp.Ncomp()) {
        MSG_WARN("Contribution vector is smaller than input's number of component, excess components will not contribute");
        for (int i = contrib.size(); i < inp.Ncomp(); i++ ) contrib.push_back(false);
    }

    //compute the density of each component of inp individually
    CompFunction<D> component_density(1, inp.Ncomp());
    component_density.func_ptr->data = inp.func_ptr->data; 
    // component_density.alloc(inp.Ncomp(), true);
    copy_grid(component_density, inp);
    multiply(prec, component_density, 1.0, inp, inp, -1, false, false, true);
    
    //collect each component's density into a temporary density, which could be defined complex for the purpose of conversion.
    CompFunction<3> rho_tmp(1, 1); //one component CompFunction
    //define rho_tmp as same number field as input for allocation
    if (inp.isreal()) {
        rho_tmp.defreal();
        rho_tmp.func_ptr->data.iscomplex = 0;
    }
    if (inp.iscomplex()) {
        rho_tmp.defcomplex();
        rho_tmp.func_ptr->data.isreal = 0;
    }
    if (rho_tmp.isreal() and rho_tmp.iscomplex()) MSG_ABORT("INPUT IS BOTH REAL AND COMPLEX; ERROR");
    rho_tmp.alloc(1, true);
    for (int i = 0; i < inp.Ncomp(); i++) {
        if (inp.isreal() and std::abs(component_density.func_ptr->data.c1[i].imag())<MachinePrec){
            if (contrib[i]) rho_tmp.CompD[0]->add_inplace((component_density.func_ptr->data.c1[i]).real(), *component_density.CompD[i]);
        } else {
            if (contrib[i]) {
                if (inp.isreal()) {
                    //this loop will at most run once per call, because we change inp.isreal()==true to ==false
                    for (int comp=0; comp < inp.Ncomp(); comp++){
                        rho_tmp.CompC[comp] = rho_tmp.CompD[comp]->CopyTreeToComplex();
                        delete rho_tmp.CompD[comp];
                        rho_tmp.CompD[comp] = nullptr; 
                    }
                    rho_tmp.defcomplex();
                }
                rho_tmp.CompC[0]->add_inplace(component_density.func_ptr->data.c1[i], *component_density.CompC[i]);
            }
        }
    }
    //making sure that out is real and clear the rest
    if (rho_tmp.iscomplex()) {
        // copy onto out's real component
        auto *realPart = rho_tmp.CompC[0]->Real(); //due to Real() returning a raw pointer, we have to create a temp variable to avoid a memory leak
        realPart->deep_copy(out.CompD[0]);
        delete realPart;
        delete rho_tmp.CompC[0];
        rho_tmp.CompC[0] = nullptr;
        out.defreal();
    } else {
        //deep copy rho_tmp's value into out
        rho_tmp.CompD[0]->deep_copy(out.CompD[0]);
        delete rho_tmp.CompD[0];
        rho_tmp.CompD[0] = nullptr;
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
        if (!out_allocated) out.alloc(out.Ncomp(), true);
        return;
    }
    if (!out_allocated) out.alloc(inp_a.Ncomp(), true);
    bool inp_aisReal = inp_a.isreal();
    bool inp_bisReal = inp_b.isreal();
    for (int comp = 0; comp < inp_a.Ncomp(); comp++) {
        if (conjugate){
            out.func_ptr->data.c1[comp] = std::conj(inp_a.func_ptr->data.c1[comp]) * inp_b.func_ptr->data.c1[comp];
        }else {
            out.func_ptr->data.c1[comp] = inp_a.func_ptr->data.c1[comp] * inp_b.func_ptr->data.c1[comp]; // we could put this is coef if everything is real?
        }
        if (inp_aisReal and inp_bisReal) {
            // if (!out_allocated) out.alloc(out.Ncomp()); //moved out of the loop
            if (need_to_multiply) {
                if (prec < 0.0) {
                    // Union grid
                    build_grid(*out.CompD[comp], *inp_a.CompD[comp]);
                    build_grid(*out.CompD[comp], *inp_b.CompD[comp]);
                    inp_a.CompD[comp]->calcSquareNorm();
                    inp_b.CompD[comp]->calcSquareNorm();
                    mrcpp::multiply(prec, *out.CompD[comp], coef, *inp_a.CompD[comp], *inp_b.CompD[comp], 0, false, false, conjugate);
                    out.CompD[comp]->calcSquareNorm();
                } else {
                    // Adaptive grid
                    inp_a.CompD[comp]->calcSquareNorm();
                    inp_b.CompD[comp]->calcSquareNorm();
                    mrcpp::multiply(prec, *out.CompD[comp], coef, *inp_a.CompD[comp], *inp_b.CompD[comp], maxIter, absPrec, useMaxNorms, conjugate);
                    out.CompD[comp]->calcSquareNorm();
                }
            }
        } else {
            // At least one of the inputs is complex
            // if one of the input is real, we simply make a new complex copy of it
            // bool inp_aisReal = inp_a.isreal();
            // bool inp_bisReal = inp_b.isreal();
            // Here, we keep the real trees of the inputs, as well as creating the complex ones, to avoid copying the trees back after multiplication. We restore the original state of the inputs after multiplication.
            if (inp_aisReal) {
                //NOTE! No need to copy every component right now as we will only use the current component
                inp_a.CompC[comp] = inp_a.CompD[comp]->CopyTreeToComplex();
                inp_a.func_ptr->iscomplex = true; //will be reverted later
                inp_a.func_ptr->isreal = false;
            }
            if (inp_bisReal) {
                inp_b.CompC[comp] = inp_b.CompD[comp]->CopyTreeToComplex();
                inp_b.func_ptr->iscomplex = true; //will be reverted later
                inp_b.func_ptr->isreal = false;
            }
            ComplexDouble coef = 1.0;
            if (need_to_multiply) {
                if (prec < 0.0) {
                    // Union grid
                    out.func_ptr->iscomplex = 1;
                    out.func_ptr->isreal = 0;
                    // if (!out_allocated) out.alloc(out.Ncomp()); //test debug moved out of loop
                    build_grid(*out.CompC[comp], *inp_a.CompC[comp]);
                    build_grid(*out.CompC[comp], *inp_b.CompC[comp]);
                    mrcpp::multiply(prec, *out.CompC[comp], coef, *inp_a.CompC[comp], *inp_b.CompC[comp], 0, false, false, conjugate);
                } else {
                    // Adaptive grid
                    if (out.CompD[comp] != nullptr) { // NB: func_ptr has alreadybeen overwritten!
                        if (out.CompD[comp]->getNNodes() > 0) {
                            out.CompC[comp] = out.CompD[comp]->CopyTreeToComplex();
                            out.func_ptr->iscomplex = 1;
                            out.func_ptr->isreal = 0;
                            delete out.CompD[comp];
                            out.CompD[comp] = nullptr;
                        } else {
                            out.func_ptr->iscomplex = 1;
                            out.func_ptr->isreal = 0;
                            // out.alloc(out.Ncomp()); //test debug
                        }
                    } else {
                        out.func_ptr->iscomplex = 1;
                        out.func_ptr->isreal = 0;
                    }
                    mrcpp::multiply(prec, *out.CompC[comp], coef, *inp_a.CompC[comp], *inp_b.CompC[comp], maxIter, absPrec, useMaxNorms, conjugate);
                }
            }
        }
    }
    // restore original tree by deleting the temporary complex tree. The real tree still exists, but is not used in the multiplication.
    if (inp_aisReal) {
        for (int comp = 0; comp < inp_a.Ncomp(); comp++) {
            delete inp_a.CompC[comp];
            inp_a.CompC[comp] = nullptr;
        }
        inp_a.func_ptr->iscomplex = false;
        inp_a.func_ptr->isreal = true;
    }
    if (inp_bisReal) {
        for (int comp = 0; comp < inp_a.Ncomp(); comp++) {
            delete inp_b.CompC[comp];
            inp_b.CompC[comp] = nullptr;
        }
        inp_b.func_ptr->iscomplex = false;
        inp_b.func_ptr->isreal = true;
    }
    mpi::share_function(out, 0, 9911, mpi::comm_share);
}

/** @brief out = inp_a * f
 *
 *  Only one component is multiplied
 */
template <int D> void multiply(CompFunction<D> &out, CompFunction<D> &inp_a, RepresentableFunction<D, double> &f, double prec, int nrefine, bool conjugate) {
    MSG_ABORT("Not implemented");
    if (inp_a.Ncomp() > 1) MSG_ABORT("Not implemented");
    if (inp_a.isreal() != 1) MSG_ABORT("Not implemented");
    if (conjugate) MSG_ABORT("Not implemented");
    CompFunctionVector CompVec; // Should use vector<CompFunction>?
    CompVec.push_back(inp_a);
    CompFunctionVector CompVecOut;
    // CompVecOut = multiply(CompVec, f, prec, nullptr, nrefine, true);
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
    // CompVecOut = multiply(CompVec, f, prec, nrefine, true); //why was this commented out? 
    out = CompVecOut[0];
}

/** @brief out = inp_a * f
 *
 */
template <int D> void multiply(CompFunction<D> &out, FunctionTree<D, double> &inp_a, RepresentableFunction<D, double> &f, double prec, int nrefine, bool conjugate) {
    CompFunction<D> func_a;
    func_a.func_ptr->isreal = 1;
    func_a.func_ptr->iscomplex = 0;
    // func_a.alloc(1);
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

/** @brief Multiplies each component of inp_a with the functionTree inp_b
 * @param out : Output compfunction
 * @param inp_a : input compfunction (e.g. a spinor)
 * @param inp_b : input FunctionTree (a scalar function)
 * @param conjugate 
 */
template <int D> void multiply(CompFunction<D> &out, CompFunction<D> inp_a, FunctionTree<D, double> &inp_b, double prec, bool absPrec, bool useMaxNorms, bool conjugate){
    if (inp_a.func_ptr->conj) conjugate = (not conjugate);
    bool need_to_multiply = not(out.isShared()) or mpi::share_master();
    bool out_allocated = true;
    if (out.Ncomp() == 0) out_allocated = false;
    bool share = out.isShared();
    out.func_ptr->data = inp_a.func_ptr->data;
    out.func_ptr->data.shared = share; // we don't inherit the shareness
    out.func_ptr->conj = false;        // we don't inherit conjugaison
    if (inp_a.getNNodes() == 0 or inp_b.getNNodes() == 0) {
        if (!out_allocated) out.alloc(inp_a.Ncomp());
        return;
    }
    double coef = 1.0;

    if (!out_allocated) out.alloc(inp_a.Ncomp()); 
    for (int comp = 0; comp < inp_a.Ncomp(); comp++) {
        if (conjugate) {
            out.func_ptr->data.c1[comp] = std::conj(inp_a.func_ptr->data.c1[comp]);
        } else{
            out.func_ptr->data.c1[comp] = inp_a.func_ptr->data.c1[comp]; // we could put this is coef if everything is real?
        }
        if (inp_a.isreal()) {
            if (need_to_multiply) {
                // if (!out_allocated) out.alloc(inp_a.Ncomp());//moved outside the loop
                if (prec < 0.0) {
                    // Union grid
                    build_grid(*out.CompD[comp], *inp_a.CompD[comp]);
                    build_grid(*out.CompD[comp], inp_b);
                    mrcpp::multiply(prec, *out.CompD[comp], coef, *inp_a.CompD[comp], inp_b, 0, false, false, conjugate); //NOTE: maybe maxIter should be -1 instead of 0?
                } else {
                    // Adaptive grid
                    mrcpp::multiply(prec, *out.CompD[comp], coef, *inp_a.CompD[comp], inp_b, -1, absPrec, useMaxNorms, conjugate);
                }
            }
        } else {
            // inp_a is complex
            // therefore we need to create a complex copy of inp_b
            FunctionTree<D, ComplexDouble> *pointer_to_inp_b_comp = nullptr; //might not be the 
            pointer_to_inp_b_comp = inp_b.CopyTreeToComplex();

            ComplexDouble coef = 1.0;
            if (need_to_multiply) {
                if (prec < 0.0) {
                    // Union grid 
                    out.func_ptr->iscomplex = 1;
                    out.func_ptr->isreal = 0;
                    delete out.CompD[comp];
                    // out.CompD[comp] = nullptr; //this should be added eventually but it makes the mrchem behave weirdly even if it never gets called
                    delete out.CompC[comp]; 
                    // out.CompC[comp] = new FunctionTree<D, ComplexDouble>(*defaultCompMRA<D>); //this should be added eventually but it makes the mrchem behave weirdly even if it never gets called
                    build_grid(*out.CompC[comp], *inp_a.CompC[comp]);
                    // build_grid(*out.CompC[comp], inp_b);
                    mrcpp::multiply(prec, *out.CompC[comp], coef, *inp_a.CompC[comp], *pointer_to_inp_b_comp, 0, false, false, conjugate);
                } else { // note that this assumes Ncomp=1
                    // Adaptive grid
                    if (out.CompD[comp] != nullptr) { // NB: func_ptr has alreadybeen overwritten!
                        if (out.CompD[comp]->getNNodes() > 0) {
                            out.CompC[comp] = out.CompD[comp]->CopyTreeToComplex();
                            out.func_ptr->iscomplex = 1;
                            out.func_ptr->isreal = 0;
                            delete out.CompD[comp];
                            out.CompD[comp] = nullptr;
                        } else {
                            out.func_ptr->iscomplex = 1;
                            out.func_ptr->isreal = 0;
                        }
                    } else {
                        out.func_ptr->iscomplex = 1;
                        out.func_ptr->isreal = 0;
                    }
                    mrcpp::multiply(prec, *out.CompC[comp], coef, *inp_a.CompC[comp], *pointer_to_inp_b_comp, -1, absPrec, useMaxNorms, conjugate);
                }
            }
            //free the temporary ComplexDouble copy of inp_b to avoid memory leaks
            delete pointer_to_inp_b_comp; 
            pointer_to_inp_b_comp = nullptr; //unnecessary but feels right
        }
    }
    mpi::share_function(out, 0, 9911, mpi::comm_share);
}


/*
* @brief Component-wise multiplication from the right with a complex function tree.
* @param conjugate: boolean to choose if we want to conjugate the CompFunction inp_a, not the FunctionTree.
*/
template <int D> void multiply(CompFunction<D> &out, CompFunction<D> inp_a, FunctionTree<D, ComplexDouble> &inp_b, double prec, bool absPrec, bool useMaxNorms, bool conjugate){
    if (inp_a.func_ptr->conj) conjugate = (not conjugate);
    bool need_to_multiply = not(out.isShared()) or mpi::share_master();
    bool out_allocated = true;
    if (out.Ncomp() == 0) out_allocated = false;
    bool share = out.isShared();
    out.func_ptr->data = inp_a.func_ptr->data;
    out.func_ptr->data.shared = share; // we don't inherit the shareness
    out.func_ptr->conj = false;        // we don't inherit conjugaison
    if (inp_a.getNNodes() == 0 or inp_b.getNNodes() == 0) {
        if (!out_allocated) out.alloc(inp_a.Ncomp());
        return;
    }

    double coef = 1.0;
    if (!out_allocated) out.alloc(inp_a.Ncomp(), true);
    for (int comp = 0; comp < inp_a.Ncomp(); comp++) {
        // out.func_ptr->data.c1[comp] = inp_a.func_ptr->data.c1[comp] ; // we could put this is coef if everything is real?
        if (conjugate) {
            out.func_ptr->data.c1[comp] = std::conj(inp_a.func_ptr->data.c1[comp]);
        } else{
            out.func_ptr->data.c1[comp] = inp_a.func_ptr->data.c1[comp]; // we could put this is coef if everything is real?
        }
        // Whether or not inp_a is complex, inp_b is, so we simply make a complex copy of inp_a if it is real
        // bool inp_bisReal = inp_b.isreal();
        bool inp_aisReal = inp_a.isreal(); //book keeping to restore inp_a's properties at the end of the iteration
        if (inp_aisReal) {
            inp_a.CompC[comp] = inp_a.CompD[comp]->CopyTreeToComplex();
            inp_a.func_ptr->iscomplex = true;
            inp_a.func_ptr->isreal = false;
        }
        ComplexDouble coef = 1.0;
        if (need_to_multiply) {
            if (prec < 0.0) {
                // Union grid
                out.func_ptr->iscomplex = 1;
                out.func_ptr->isreal = 0;
                delete out.CompD[comp];
                // out.CompD[comp] = nullptr; //this should be added eventually but it makes the mrchem behave weirdly even if it never gets called
                delete out.CompC[comp];  
                // out.CompC[comp] = new FunctionTree<D, ComplexDouble>(*defaultCompMRA<D>);//this should be added eventually but it makes the mrchem behave weirdly even if it never gets called
                // if (!out_allocated) out.alloc(inp_a.Ncomp());
                build_grid(*out.CompC[comp], *inp_a.CompC[comp]);
                // build_grid(*out.CompC[comp], inp_b);
                mrcpp::multiply(prec, *out.CompC[comp], coef, *inp_a.CompC[comp], inp_b, 0, false, false, conjugate);
            } else { // note that this assumes Ncomp=1 -- still true? 
                // Adaptive grid
                if (out.CompD[comp] != nullptr) { // NB: func_ptr has alreadybeen overwritten!
                    if (out.CompD[comp]->getNNodes() > 0) {
                        out.CompC[comp] = out.CompD[comp]->CopyTreeToComplex();
                        out.func_ptr->iscomplex = 1;
                        out.func_ptr->isreal = 0;
                        delete out.CompD[comp];
                        out.CompD[comp] = nullptr;
                    } else {
                        out.func_ptr->iscomplex = 1;
                        out.func_ptr->isreal = 0;
                        // out.alloc(inp_a.Ncomp());
                    }
                } else {
                    out.func_ptr->iscomplex = 1;
                    out.func_ptr->isreal = 0;
                    // if (!out_allocated) out.alloc(inp_a.Ncomp());
                }
                mrcpp::multiply(prec, *out.CompC[comp], coef, *inp_a.CompC[comp], inp_b, -1, absPrec, useMaxNorms, conjugate);
            }
        }
        // restore original tree
        if (inp_aisReal) {
            delete inp_a.CompC[comp];
            inp_a.CompC[comp] = nullptr;
            inp_a.func_ptr->iscomplex = false;
            inp_a.func_ptr->isreal = true;
        }
        // }
    }
    mpi::share_function(out, 0, 9911, mpi::comm_share);
}

/** @brief Compute <bra|ket> = \int_Domain \bra^\dagger(r) * \ket(r) dr.
 *
 *  Sum of component dots.
 *  Notice that the <bra| position is complex conjugated in the tree multiplication.
 *
 *  NOTE: if the number of components do not match, the minimum of the two is taken. 
 *  Doesn't throw an error or warning (despite making no mathematical sense), 
 *  because it happens a lot when using MPI.
 */
// template <int D> ComplexDouble dot(CompFunction<D> bra, CompFunction<D> ket) {
template <int D> ComplexDouble dot(const CompFunction<D> &bra, const CompFunction<D> &ket) {
    // if (bra.Ncomp() != ket.Ncomp()) MSG_WARN("Mismatched number of components: Ncomp_bra="<< bra.Ncomp() << ", Ncomp_ket="<< ket.Ncomp() << ", taking the minimum between the two");
    ComplexDouble dotprodtot = 0.0;
    for (int comp = 0; comp < std::min(bra.Ncomp(), ket.Ncomp()); comp++) {
        ComplexDouble dotprod = 0.0;
        //Computing the dot product of the current components for each
        //case of bra/ket being real or complex
        if (bra.isreal() and ket.isreal()) {
            dotprod += mrcpp::dot(*bra.CompD[comp], *ket.CompD[comp]);
        } else if (bra.isreal() and ket.iscomplex()) {
            dotprod += mrcpp::dot(*bra.CompD[comp], *ket.CompC[comp]);
        } else if (bra.iscomplex() and ket.isreal()) {
            dotprod += mrcpp::dot(*bra.CompC[comp], *ket.CompD[comp]);
        } else { //both complex
            dotprod += mrcpp::dot(*bra.CompC[comp], *ket.CompC[comp]);
        }
        //Multiplying the dot product of the current components by the c1 coefficients
        dotprod *= std::conj(bra.func_ptr->data.c1[comp]) * ket.func_ptr->data.c1[comp];
        //Adding the result to the total dot product
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

/* @brief Projects real f onto out 
    * @param out Output CompFunction
    * @param f Input function to project
    * @param prec Precision for the projection
    * @param comp Component index to project onto
*/
void project(CompFunction<3> &out, std::function<double(const Coord<3> &r)> f, double prec, int comp) {
    //mpi shenanigans
    bool need_to_project = not(out.isShared()) or mpi::share_master();
    out.func_ptr->isreal = 1;
    out.func_ptr->iscomplex = 0;
    // allocating a component if compFunction is empty
    if (out.Ncomp() < 1) out.alloc(1);
    // This variable ensures that scalar function do not allocate additional components
    bool need_to_allocate = true;
    if (out.Ncomp() < 2) {need_to_allocate = false;}

    // lambda function when we need to initialise a component to zero
    std::function<double(const Coord<3>&)> fzero = [](const Coord<3> &r) -> double{ return 0.0; };

    // allocating and projecting the component(s)
    for (int i = 0; i < out.Ncomp(); i++) {
        if (i == comp) {
            out.alloc_comp(i);
            mrcpp::project<3>(prec, *out.CompD[i], f);
        } else if (need_to_allocate) {
            // out.CompD[i]->setZero();
            out.alloc_comp(i, true);
            // mrcpp::project<3>(prec, *out.CompD[i], fzero);
        }
    };
    //mpi 
    mpi::share_function(out, 0, 132231, mpi::comm_share);
}

void project_real(CompFunction<3> &out, std::function<double(const Coord<3> &r)> f, double prec) {
    bool need_to_project = not(out.isShared()) or mpi::share_master();
    out.func_ptr->isreal = 1;
    out.func_ptr->iscomplex = 0;
    if (out.Ncomp() < 1) out.alloc(1);
    for (int comp=0; comp < out.Ncomp(); comp++) if(out.CompD[comp]==nullptr) out.alloc_comp(comp, true);
    if (need_to_project) mrcpp::project<3>(prec, *out.CompD[0], f);
    mpi::share_function(out, 0, 123123, mpi::comm_share);
}

// template <int D, typename T>
/* @brief Projects complex f onto out 
    * @param out Output CompFunction
    * @param f Input function to project
    * @param prec Precision for the projection
    * @param comp Component index to project onto
*/
void project(CompFunction<3> &out, std::function<ComplexDouble(const Coord<3> &r)> f, double prec, int comp) {
    bool need_to_project = not(out.isShared()) or mpi::share_master();
    out.func_ptr->isreal = 0;
    out.func_ptr->iscomplex = 1;

    if (out.Ncomp() < 1) out.alloc(1); // probably should just return an error 
    // This variable ensures that scalar function do not allocate additional components and initialise them to zero
    bool need_to_allocate = true;
    if (out.Ncomp() < 2) {need_to_allocate = false;}

    // lambda function when we need to initialise a component to zero
    std::function<ComplexDouble(const Coord<3>&)> fzero = [](const Coord<3> &r) -> ComplexDouble { return ComplexDouble(0.0, 0.0); };

    // allocating and projecting the component(s)
    for (int i = 0; i < out.Ncomp(); i++) {
        if (i == comp) {
            out.alloc_comp(i);
            mrcpp::project<3>(prec, *out.CompC[i], f);
        } else if (need_to_allocate) { //This may be completely useless and might be replaceable by setZero, but the code seems to be working and I don't dare break it now
            // out.CompC[i]->setZero();
            out.alloc_comp(i, true);
            // mrcpp::project<3>(prec, *out.CompC[i], fzero);
        }
    };
    mpi::share_function(out, 0, 123123, mpi::comm_share);
}

void project_cplx(CompFunction<3> &out, std::function<ComplexDouble(const Coord<3> &r)> f, double prec) {
    bool need_to_project = not(out.isShared()) or mpi::share_master();
    out.func_ptr->isreal = 0;
    out.func_ptr->iscomplex = 1;
    if (out.Ncomp() < 1) out.alloc(1);
    for (int comp=0; comp < out.Ncomp(); comp++) if(out.CompC[comp]==nullptr) out.alloc_comp(comp, true);
    if (need_to_project) mrcpp::project<3>(prec, *out.CompC[0], f);
    mpi::share_function(out, 0, 123123, mpi::comm_share);
}


/** @brief Project a RepresentableFunction onto a real-valued CompFunction
 * @param out Output CompFunction
 * @param f Input RepresentableFunction 
 * @param prec Precision for the projection
 * @param nComp Number of components to project onto
 */
template <int D> void project(CompFunction<D> &out, RepresentableFunction<D, double> &f, double prec, int nComp) {
    //mpi ownership check
    bool need_to_project = not(out.isShared()) or mpi::share_master();
    //real case
    out.defreal();

    //free + reallocate out. Its content (if any) was going to get trashed anyway
    if (out.Ncomp() < 1) {
        // allocating a component if compFunction is empty or if more components than are present must be projected
        out.alloc(nComp, true);
        out.func_ptr->data.Ncomp = nComp;
    } else { // simply freeing the component(s)
        out.alloc(out.Ncomp(), true);
    }
    
    // allocating and projecting the component(s)
    if (need_to_project){
        for(int c=0; c<nComp; c++) { // (;ºoº;) [c++] ¬-("o")   
            build_grid(*out.CompD[c], f);
            mrcpp::project<D, double>(prec, *out.CompD[c], f);
        } 
    }

    //mpi
    mpi::share_function(out, 0, 132231, mpi::comm_share);
}

/** @brief Project a complex RepresentableFunction onto a complex-defined CompFunction
 * @param out Output CompFunction
 * @param f Input RepresentableFunction 
 * @param prec Precision for the projection
 * @param nComp Number of components to project onto
 */
template <int D> void project(CompFunction<D> &out, RepresentableFunction<D, ComplexDouble> &f, double prec, int nComp) {
    //mpi ownership check
    bool need_to_project = not(out.isShared()) or mpi::share_master();
    //complex case
    out.defcomplex();

    //free + reallocate out. Its content (if any) was going to get trashed anyway
    if (out.Ncomp() < 1) {
        // allocating a component if compFunction is empty or if more components than are present must be projected
        out.alloc(nComp, true);
        out.func_ptr->data.Ncomp = nComp;
    } else { // simply freeing the component(s)
        out.alloc(out.Ncomp(), true);
    }
    
    // allocating and projecting the component(s)
    if (need_to_project){
        for(int c=0; c<nComp; c++) { // (;ºoº;) [c++] ¬-("o")   
            build_grid(*out.CompC[c], f);
            mrcpp::project<D, ComplexDouble>(prec, *out.CompC[c], f);
        } 
    }

    //mpi
    mpi::share_function(out, 0, 132231, mpi::comm_share);
}



// =============================================================================
// =========================CompFunctionVector==================================
// =============================================================================

CompFunctionVector::CompFunctionVector(int N)
        : std::vector<CompFunction<3>>(N) {
    for (int i = 0; i < N; i++) (*this)[i].func_ptr->rank = i;
    vecMRA = defaultCompMRA<3>;
}
void CompFunctionVector::distribute() {
    for (int i = 0; i < this->size(); i++) (*this)[i].func_ptr->rank = i;
}

// CompFunction<3> CompFunctionVector::operator[](int i) const {
//     if (i < 0 || i >= this->size()) {
//         throw std::out_of_range("Index out of range in CompFunctionVector");
//     }
//     return this->at(i); 
// }

// void project(CompFunctionVector &out, std::function<double(const Coord<3> &r)> f, int index double prec) {
//     bool need_to_project = not(out.isShared()) or mpi::share_master();
//     for (int i = 0; i < out.size(); i++) {
//         out[i].func_ptr->isreal = 1;
//         out[i].func_ptr->iscomplex = 0;
//         if (out[i].Ncomp() < 1) out[i].alloc(1);
//     }
//     if (need_to_project) mrcpp::project<3>(prec, out, f);
//     mpi::share_function(out, 0, 123123, mpi::comm_share);
// }

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
    for (int i = 0; i < N; i++) {
        if (!mrcpp::mpi::my_func(i)) continue;
        if (Phi[i].isreal()){
            for (int j = 0; j < Phi[i].Ncomp(); j++) {
                Phi[i].CompC[j] = Phi[i].CompD[j]->CopyTreeToComplex();
                delete Phi[i].CompD[j];
                Phi[i].CompD[j] = nullptr;
            }
        }
        Phi[i].func_ptr->isreal = 0;
        Phi[i].func_ptr->iscomplex = 1;
    }
    int M = Psi.size();
    for (int i = 0; i < M; i++) {
        if (!mrcpp::mpi::my_func(i)) continue;
        for (int j = 0; j < 4; j++) {
            if (Psi[i].CompD[j] != nullptr){
                delete Psi[i].CompD[j];
                Psi[i].CompD[j] = nullptr;
            }
        }
        Psi[i].func_ptr->isreal = 0;
        Psi[i].func_ptr->iscomplex = 1;
    }

    for (int i = 0; i < N; i++) {
        if (!mrcpp::mpi::my_func(i)) continue;
        if (Phi[i].func_ptr->conj) MSG_ABORT("Conjugaison not implemented for rotations");
    }
    if (U.rows() < N) MSG_ABORT("Incompatible number of rows for U matrix");
    if (U.cols() < M) MSG_ABORT("Incompatible number of columns for U matrix");

    // Compute the number of components from a function this MPI worker owns, to avoid issues
    int Ncomponents = 1; 
    for (int i = 0; i < N; i++) {
        if (!mrcpp::mpi::my_func(i)) continue;
        Ncomponents = std::max(Ncomponents, Phi[i].Ncomp());
    }
    // start rotating
    for (int q = 0; q < Ncomponents; q++) {    
        // 1) make union tree without coefficients. Note that the ref tree is always real (in fact it has no coeff)
        FunctionTree<3> refTree(*Phi.vecMRA);
        mpi::allreduce_Tree_noCoeff(refTree, Phi, mpi::comm_wrk, q);

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
            if (!mpi::my_func(i)) continue;
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
                Phi[j].complex(q).makeCoeffVector(coeffVec[j], indexVec[j], parindexVec, scalefac, max_ix, refTree);
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
            save_nodes(Phi, refTree, nodesPhi, -1, q); //only current component at once
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
                // if (node2orbVec[node_ix].size() <= 0) continue;
                auto it = node2orbVec.find(node_ix); //test debug test multithreading
                if (it == node2orbVec.end() || it->second.empty()) continue; //test debug test multithreading
                csize = sizecoeffW;
                if (parindexVec_ref[n] < 0) csize = sizecoeff; // for root nodes we include scaling coeff

                int shift = sizecoeff - sizecoeffW; // to copy only wavelet part
                if (parindexVec_ref[n] < 0) shift = 0;
                ComplexMatrix coeffBlock(csize, node2orbVec[node_ix].size());
                // for (int j : node2orbVec[node_ix]) { // loop over indices of the orbitals using this node
                for (int j : it->second) { // loop over indices of the orbitals using this node
                    // int orb_node_ix = orb2node[j][node_ix];
                    int orb_node_ix = orb2node[j].at(node_ix);
                    for (int k = 0; k < csize; k++) coeffBlock(k, orbjVec.size()) = coeffVec[j][orb_node_ix][k + shift];
                    orbjVec.push_back(j);
                }

                // 4b) make a list of rotated orbitals needed for this node
                // OMP must wait until parent is ready
                // while (parindexVec_ref[n] >= 0 and nodeReady[ix2coef_ref[parindexVec_ref[n]]] == 0) {
                while (parindexVec_ref[n] >= 0 and nodeReady[ix2coef_ref.at(parindexVec_ref[n])] == 0) {
    #pragma omp flush
                };

                std::vector<int> orbiVec;
                for (int i = 0; i < M; i++) {                                                                        // loop over all rotated orbitals
                    // if (parindexVec_ref[n] >= 0 and split_serial(i, ix2coef_ref[parindexVec_ref[n]]) == 0) continue; // parent node has too small wavelets
                    if (parindexVec_ref[n] >= 0 and split_serial(i, ix2coef_ref.at(parindexVec_ref[n])) == 0) continue; // parent node has too small wavelets
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
                    // if (parindexVec_ref[n] == -1 or split_serial(orbiVec[i], ix2coef_ref[parindexVec_ref[n]])) {
                    if (parindexVec_ref[n] == -1 or split_serial(orbiVec[i], ix2coef_ref.at(parindexVec_ref[n]))) {
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
    #pragma omp critical //preventive construct to avoid race conditions, because it appears in the if statement so might as well copy it here? 
                        ix2coef[orbiVec[i]][node_ix] = max_n + 1; // should not be used
                        split_serial(orbiVec[i], n) = 0;          // do not split if parent does not need to be split
                    }
                }
    #pragma omp flush //preventive construct to avoid race conditions (useless?)
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
                // Psi[j].alloc(1); // All data is stored in coeffpVec[j]
                Psi[j].alloc_comp(q);
                Psi[j].complex(q).clear();
                Psi[j].complex(q).makeTreefromCoeff(refTree, coeffpVec[j], ix2coef[j], prec);
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

                // Psi[j].alloc(1);
                Psi[j].alloc_comp(q); //could be a problem to allocate in parallel, but it should be okay since there is no OMP here?
                Psi[j].complex(q).clear();
                Psi[j].complex(q).makeTreefromCoeff(refTree, coeffpVec, ix2coef, prec);

                for (ComplexDouble *p : pointerstodelete) delete[] p;
                pointerstodelete.clear();
            }
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
    bool iscomplex = false;
    int N = Phi.size();
    int M = Psi.size();

    // The principle of this routine is that nodes are rotated one by one using matrix multiplication.
    // The routine does avoid when possible to move data, but uses pointers and indices manipulation.
    // MPI version does not use OMP yet, Serial version uses OMP
    // size of input is N, size of output is M
    if (U.rows() < N) MSG_ABORT("Incompatible number of rows for U matrix");
    if (U.cols() < M) MSG_ABORT("Incompatible number of columns for U matrix");

    // Rescaling the trees with their prefactors c1 (which might make them complex, hence why we do it here)
    // we first syncronize all the number of components
    // we assume that at least one orbital is owned by this MPI (TODO: allreduce)

    for (int i=0; i < N; i++){
        if (!mrcpp::mpi::my_func(i)) continue;
        for (int q=0; q<Phi[i].Ncomp();q++){
            if (Phi[i].isreal() and (std::abs(Phi[i].func_ptr->data.c1[q].imag())<MachineZero)) {
                Phi[i].CompD[q]->rescale((Phi[i].func_ptr->data.c1[q]).real());
            } else {
                if (Phi[i].isreal()) {
                    for (int comp=0; comp<Phi[i].Ncomp();comp++){
                        Phi[i].CompC[comp] = Phi[i].CompD[comp]->CopyTreeToComplex();
                        delete Phi[i].CompD[comp];
                        Phi[i].CompD[comp] = nullptr; 
                    }
                    Phi[i].defcomplex();
                }
                //Now that the trees are converted, rescaling
                Phi[i].CompC[q]->rescale(Phi[i].func_ptr->data.c1[q]);
                //Resetting the prefactor to 1, as it is now held in the tree itself
                Phi[i].func_ptr->data.c1[q] = {1.0,0.0};
            }
        }
    }

    //Handling complex case
    for (int i=0; i < N; i++) {
        //checking if the orbitals are complex
        if (not mrcpp::mpi::my_func(i)) continue;
        if (Phi[i].iscomplex()) iscomplex=true;
    }
    //checking if the rotation matrix is complex-valued
    iscomplex += (U.imag().cwiseAbs().maxCoeff() > mrcpp::MachineZero); //Cursed code brought to you by Niklas
    //sync across mpi ranks
    iscomplex = mrcpp::mpi::allreduce_max(iscomplex ? 1 : 0, mrcpp::mpi::comm_wrk) > 0;

    if (iscomplex) {
        rotate_cplx(Phi, U, Psi, prec);
        return;
    }

    // Computing the number of components of spinors in an MPI safe way
    int Ncomponents = 1;
    for (int i = 0; i < N; i++) {
        if (!mrcpp::mpi::my_func(i)) continue;
        Ncomponents = std::max(Ncomponents, Phi[i].Ncomp());
    }
    for (int q = 0; q < Ncomponents; q++) {
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
            if (!mrcpp::mpi::my_func(i)) continue;
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
                Phi[j].real(q).makeCoeffVector(coeffVec[j], indexVec[j], parindexVec, scalefac, max_ix, refTree); //adapt, similar to calc_overlap_matrix
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
            save_nodes(Phi, refTree, nodesPhi, -1, q);
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
                // if (node2orbVec[node_ix].size() <= 0) continue;
                // test debug test start
                auto it = node2orbVec.find(node_ix);
                if (it == node2orbVec.end() || it->second.empty()) continue;
                // test debug test end
                csize = sizecoeffW;
                if (parindexVec_ref[n] < 0) csize = sizecoeff; // for root nodes we include scaling coeff

                int shift = sizecoeff - sizecoeffW; // to copy only wavelet part
                if (parindexVec_ref[n] < 0) shift = 0;
                DoubleMatrix coeffBlock(csize, node2orbVec[node_ix].size());
                // for (int j : node2orbVec[node_ix]) { // loop over indices of the orbitals using this node
                for (int j : it->second) { // test debug test loop over indices of the orbitals using this node
                    // int orb_node_ix = orb2node[j][node_ix];
                    int orb_node_ix = orb2node[j].at(node_ix); // test debug test
                    for (int k = 0; k < csize; k++) coeffBlock(k, orbjVec.size()) = coeffVec[j][orb_node_ix][k + shift];
                    orbjVec.push_back(j);
                }

                // 4b) make a list of rotated orbitals needed for this node
                // OMP must wait until parent is ready
                // while (parindexVec_ref[n] >= 0 and nodeReady[ix2coef_ref[parindexVec_ref[n]]] == 0) {
                while (parindexVec_ref[n] >= 0 and nodeReady[ix2coef_ref.at(parindexVec_ref[n])] == 0) {
        #pragma omp flush
                };
                std::vector<int> orbiVec;
                for (int i = 0; i < M; i++) {                                                                        // loop over all rotated orbitals
                    // if (parindexVec_ref[n] >= 0 and split_serial(i, ix2coef_ref[parindexVec_ref[n]]) == 0) continue; // parent node has too small wavelets
                    if (parindexVec_ref[n] >= 0 and split_serial(i, ix2coef_ref.at(parindexVec_ref[n])) == 0) continue; // parent node has too small wavelets
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
                    // if (parindexVec_ref[n] == -1 or split_serial(orbiVec[i], ix2coef_ref[parindexVec_ref[n]])) {
                    if (parindexVec_ref[n] == -1 or split_serial(orbiVec[i], ix2coef_ref.at(parindexVec_ref[n]))) {
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
                        #pragma omp critical //preventive construct to avoid race conditions -- Maybe useless
                        ix2coef[orbiVec[i]][node_ix] = max_n + 1; // should not be used
                        split_serial(orbiVec[i], n) = 0;          // do not split if parent does not need to be split
                    }
                }
        #pragma omp flush //preventive construct to avoid race conditions (useless?)
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
        // std::cout << "mrcpp::CompFunction::rotate pre-reconstruction component " << q << " with Phi (1st argument) norm " << Phi[0].norm() << " with Psi norm " << Psi[0].norm() << std::endl;


        // 5) reconstruct trees using rotated nodes.

        // only serial case can use OMP, because MPI cannot be used by threads
        if (serial) {
            // OMP parallelized, but does not scale well, because the total memory bandwidth is a bottleneck. (the main
            // operation is writing the coefficient into the tree)
        
            //debug
            // CompFunction test(*defaultCompMRA<3>, 1); 
            // test.defreal();
            // test.real(0).makeTreefromCoeff(refTree, coeffpVec[0], ix2coef[0], prec);
            // std::cout << "Test reconstruction norm " << test.norm() << std::endl;
        #pragma omp parallel for schedule(static)
        // #pragma omp critical // test debug test
            for (int j = 0; j < M; j++) {
                if (coeffpVec[j].size() == 0) continue;
                // Psi[j].alloc(1);
                Psi[j].alloc_comp(q, false); //possible issue (multithreading) but maybe safe since no MPI here? 
                Psi[j].real(q).clear();
                Psi[j].real(q).makeTreefromCoeff(refTree, coeffpVec[j], ix2coef[j], prec);
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
                // Psi[j].alloc(1);
                Psi[j].alloc_comp(q, false);
                Psi[j].real(q).makeTreefromCoeff(refTree, coeffpVec, ix2coef, prec);

                for (double *p : pointerstodelete) delete[] p;
                pointerstodelete.clear();
            }
        }
    }

}

void rotate(CompFunctionVector &Phi, const ComplexMatrix &U, double prec) {
    //create a deep_copy of the input, because the complex instance of rotate 
    //erase the output's trees, which creates a seg fault when it is called
    //through this function.
    CompFunctionVector Psi(Phi.size());
    for (int i = 0; i < Phi.size(); i++) {
        Psi[i] = Phi[i].paramCopy(); //shallow copy
        deep_copy(Psi[i], Phi[i]); //make sure it is different
    }
    rotate(Psi, U, Phi, prec); 
    //Note on MPI calls here for future reference (written by Claude):
    // Full-size, index-preserving copy: deep_copy is a cheap no-op on
    // non-owned (placeholder) slots, so this needs no MPI communication.
    return;
}

/** @brief Save all nodes in bank; identify them using serialIx from refTree
 * shift is a shift applied in the id
 * @param comp: target component of the CompFunctions. MPI ranks CompFunctions in their entirety, so all components belong to the same rank.
 */
void save_nodes(CompFunctionVector &Phi, FunctionTree<3> &refTree, BankAccount &account, int sizes, int comp) {
    int sizecoeff = (1 << refTree.getDim()) * refTree.getKp1_d();
    int sizecoeffW = ((1 << refTree.getDim()) - 1) * refTree.getKp1_d();
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
            Phi[j].real(comp).makeCoeffVector(coeffVec, indexVec, parindexVec, scalefac, max_ix, refTree);
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
            Phi[j].complex(comp).makeCoeffVector(coeffVec_cplx, indexVec, parindexVec, scalefac, max_ix, refTree);
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
            if (node2orbVec.count(node_ix) == 0) continue;

            // 3a) make values for f at this node
            // 3a1) get coordinates of quadrature points for this node
            Eigen::MatrixXd pts; // Eigen::Zero(D, nCoefs);
            std::vector<double> fval(nCoefs);
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
                    Fnode->attachCoefs(fval.data()); // note that each thread has its own copy
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
        TaskManager tasks(max_n);
        for (int nn = 0; nn < max_n; nn++) {
            int n = tasks.next_task();
            if (n < 0) break;
            MWNode<D> node(*(refNodes[n]), false);
            // 3a) make values for f
            // 3a1) get coordinates of quadrature points for this node
            Eigen::MatrixXd pts;           // Eigen::Zero(D, nCoefs);
            node.getExpandedChildPts(pts); // TODO: use getPrimitiveChildPts (less cache).
            std::vector<double> fval(nCoefs);
            Coord<D> r;
            MWNode<D> Fnode(*(refNodes[n]), false);
            if (Func == nullptr) {
                for (int j = 0; j < nCoefs; j++) {
                    for (int d = 0; d < D; d++) r[d] = pts(d, j); //*scaling_factor[d]?
                    fval[j] = f.evalf(r);
                }
            } else {
                int nIdx = Func->real().getIx(node.getNodeIndex());
                if (nIdx < 0) {
                    // use the function f instead of Func
                    for (int j = 0; j < nCoefs; j++) {
                        for (int d = 0; d < D; d++) r[d] = pts(d, j);
                        fval[j] = f.evalf(r);
                    }
                } else {
                    Func->real().getNodeCoeff(nIdx, fval.data()); // fetch coef from Bank
                    Fnode.attachCoefs(fval.data());
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
    defaultCompMRA<3> = std::make_shared<MultiResolutionAnalysis<3>>(*MRA);
}

// ComplexVector dot(CompFunctionVector &Bra, CompFunctionVector &Ket) {
//     int N = Bra.size();
//     ComplexVector result = ComplexVector::Zero(N);
//     for (int i = 0; i < N; i++) {
//         // The bra is sent to the owner of the ket
//         if (my_func(Bra[i]) != my_func(Ket[i])) { MSG_ABORT("same indices should have same ownership"); }
//         result[i] = dot(Bra[i], Ket[i]);
//         if (not mrcpp::mpi::my_func(i)) Bra[i].free();
//     }
//     mrcpp::mpi::allreduce_vector(result, mrcpp::mpi::comm_wrk);
//     return result;
// }

ComplexVector dot(CompFunctionVector &Bra, CompFunctionVector &Ket) {
    int N = Bra.size();
    ComplexVector result = ComplexVector::Zero(N);
    for (int i = 0; i < N; i++) {
        // if(!my_func(i))  result[i] = 0.0;
        // else {
        result[i] = dot(Bra[i], Ket[i]);
        // }
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

/** @brief Compute Löwdin orthonormalization matrix for a 2 component orbital vector
 *
 * @param Phi: orbitals to orthonomalize
 *
 * Computes the inverse square root of the orbital overlap matrix S^(-1/2)sk
 */
// I don't know what this would be useful for, maybe remove it?
ComplexMatrix calc_lowdin_matrix_2c(CompFunctionVector &Phi_top, CompFunctionVector &Phi_bottom) {
    MSG_ABORT("THIS FUNCTION SHOULD BE DEPRECATED - ABORTING");
    ComplexMatrix S_tilde_t = calc_overlap_matrix(Phi_top);
    ComplexMatrix S_tilde_b = calc_overlap_matrix(Phi_bottom);

    

    ComplexMatrix S_tilde = S_tilde_t + S_tilde_b; // S = S_top + S_bottom
    
    // Complex conjugate the elements of S_tilde
    for (int i = 0; i < S_tilde.rows(); i++) {
        for (int j = 0; j < S_tilde.cols(); j++) {
            S_tilde(i, j) = std::conj(S_tilde(i, j));
        }
    }


    ComplexMatrix S_m12 = math_utils::hermitian_matrix_pow(S_tilde, -0.5);
    return S_m12;
}

/** @brief Computes the overlap matrix S for complex-valued orbital vector
 *
 * NOTE: OrbitalVector is considered a ROW vector, so rotation
 *       means matrix multiplication from the right
 *
 * MPI: Rank distribution of output vector is the same as input vector
 *
 */
ComplexMatrix calc_overlap_matrix_cplx(CompFunctionVector &BraKet) {
    int N = BraKet.size();
    ComplexMatrix Stot = ComplexMatrix::Zero(N, N);
    // bool braketisreal = false;
    //simplest thing to do is to just convert every tree to complex and get rid of the real
    for (int k = 0; k < N; k++) {
        if (!mrcpp::mpi::my_func(k)) continue;
        if (BraKet[k].isreal()){
            // braketisreal = true;
            for (int comp = 0; comp< BraKet[k].Ncomp(); comp++) {
                BraKet[k].CompC[comp] = BraKet[k].CompD[comp]->CopyTreeToComplex();
                delete BraKet[k].CompD[comp];
                BraKet[k].CompD[comp] = nullptr;
            }
            BraKet[k].defcomplex();
            BraKet[k].func_ptr->data.isreal = false;
        }
    }
    // DoubleMatrix Sreal = Stot.real();
    std::shared_ptr<MultiResolutionAnalysis<3>> mra = BraKet.vecMRA;

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

    
    bool serial = mrcpp::mpi::wrk_size == 1; // flag for serial/MPI switch
    
    // Computing the number of components of the owned elements
    int Ncomponents = 1;
    for (int i = 0; i < N; i++) {
        if (!mrcpp::mpi::my_func(i)) continue;
        Ncomponents = std::max(Ncomponents, BraKet[i].Ncomp());
    }
    for (int l = 0; l < Ncomponents; l++) {
        // only used for serial case:
        std::vector<std::vector<ComplexDouble *>> coeffVec(N);
        std::map<int, std::vector<int>> node2orbVec; // for each node index, gives a vector with the indices of the orbitals using this node
        std::vector<std::map<int, int>> orb2node(N); // for a given orbital and a given node, gives the node index in
                                                    // the orbital given the node index in the reference tree
                                                    
        mrcpp::BankAccount nodesBraKet;
        ComplexMatrix S = ComplexMatrix::Zero(N, N); // Overlap matrix for this specific component, to be accumulated into Stot
        // In the serial case we store the coeff pointers in coeffVec. In the mpi case the coeff are stored in the bank
        if (serial) {
            // 2) make list of all coefficients, and their reference indices
            // for different orbitals, indexVec will give the same index for the same node in space
            std::vector<int> parindexVec; // serialIx of the parent nodes
            std::vector<int> indexVec;    // serialIx of the nodes
            for (int j = 0; j < N; j++) {
                // make vector with all coef pointers and their indices in the union grid
                BraKet[j].complex(l).makeCoeffVector(coeffVec[j], indexVec, parindexVec, scalefac, max_ix, refTree);
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
            save_nodes(BraKet, refTree, nodesBraKet, -1, l);
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
                std::map<int, std::vector<int>>::const_iterator it;
                if (serial) {
                    it = node2orbVec.find(node_ix);
                    if (it == node2orbVec.end() || it->second.empty()) continue;
                }
                if (parindexVec_ref[n] < 0)
                    csize = sizecoeff;
                else
                    csize = sizecoeffW;

                // In the serial case we copy the coeff coeffBlock. In the mpi case coeffBlock is provided by the bank
                if (serial) {
                    int shift = sizecoeff - sizecoeffW; // to copy only wavelet part
                    if (parindexVec_ref[n] < 0) shift = 0;
                    // ComplexMatrix coeffBlock(csize, node2orbVec[node_ix].size());
                    // for (int j : node2orbVec[node_ix]) { // loop over indices of the orbitals using this node
                    ComplexMatrix coeffBlock(csize, it->second.size());
                    for (int j : it->second) { // loop over indices of the orbitals using this node
                        int orb_node_ix = orb2node[j].at(node_ix); 
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
                                S(orbVec[i], orbVec[j]) += S_temp(i, j); //Partially filled with only this rank's contribution to the component overlap matrix
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
                if (i != j) S(j, i) = std::conj(S(i, j)); // ensure exact symmetry
            }
        }
        // Assumes linearity: result is sum of all nodes contributions
        mrcpp::mpi::allreduce_matrix(S, mrcpp::mpi::comm_wrk);
        // multiply by CompFunction multiplicative factor

        ComplexVector Fac = ComplexVector::Zero(N);
        for (int i = 0; i < N; i++) {
            if (!mrcpp::mpi::my_func(BraKet[i])) continue;
            Fac[i] = BraKet[i].func_ptr->data.c1[l];
        }

        mrcpp::mpi::allreduce_vector(Fac, mrcpp::mpi::comm_wrk);
        for (int i = 0; i < N; i++) {
            if (!mrcpp::mpi::my_func(BraKet[i])) continue;
            for (int j = 0; j < N; j++) { S(i, j) *= std::conj(Fac[i]) * Fac[j]; }
        }
        //Collecting the l^th component's contribution to the overlap matrix
        Stot += S;
    }

    return Stot;
}

ComplexMatrix calc_overlap_matrix(CompFunctionVector &BraKet) {
    // NB: should be spinseparated at this point!
    int N = BraKet.size();

    // Complex case determination
    bool is_complex = false;
    for (int k = 0; k < N; k++ ){
        if (!mrcpp::mpi::my_func(BraKet[k])) continue;
        if (BraKet[k].iscomplex()) is_complex = true;
    }
    is_complex = mrcpp::mpi::allreduce_max(is_complex ? 1 : 0, mrcpp::mpi::comm_wrk) > 0;
    if (is_complex) { return calc_overlap_matrix_cplx(BraKet); }

    // If this point is reached, real case
    ComplexMatrix Stot = ComplexMatrix::Zero(N, N);

    std::shared_ptr<MultiResolutionAnalysis<3>> mra = BraKet.vecMRA;

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

    
    bool serial = mrcpp::mpi::wrk_size == 1; // flag for serial/MPI switch
    
    // Computing the number of components of the owned elements
    int Ncomponents = 1;
    for (int i = 0; i < N; i++) {
        if (!mrcpp::mpi::my_func(i)) continue;
        Ncomponents = std::max(Ncomponents, BraKet[i].Ncomp());
    }
    for (int l = 0; l < Ncomponents; l++) {
        // only used for serial case:
        std::vector<std::vector<double *>> coeffVec(N);
        std::map<int, std::vector<int>> node2orbVec; // for each node index, gives a vector with the indices of the orbitals using this node
        std::vector<std::map<int, int>> orb2node(N); // for a given orbital and a given node, gives the node index in
                                                    // the orbital given the node index in the reference tree
        mrcpp::BankAccount nodesBraKet;

        ComplexMatrix S = ComplexMatrix::Zero(N, N);
        
        // In the serial case we store the coeff pointers in coeffVec. In the mpi case the coeff are stored in the bank
        if (serial) {
            // 2) make list of all coefficients, and their reference indices
            // for different orbitals, indexVec will give the same index for the same node in space
            std::vector<int> parindexVec; // serialIx of the parent nodes
            std::vector<int> indexVec;    // serialIx of the nodes
            for (int j = 0; j < N; j++) {
                // make vector with all coef pointers and their indices in the union grid
                BraKet[j].real(l).makeCoeffVector(coeffVec[j], indexVec, parindexVec, scalefac, max_ix, refTree);
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
            save_nodes(BraKet, refTree, nodesBraKet, -1, l);
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
                // if (serial and node2orbVec[node_ix].size() <= 0) continue;
                std::map<int, std::vector<int>>::const_iterator it; // test debug test multithreading
                if (serial){ //test debug test multithreading
                    it = node2orbVec.find(node_ix); //test debug test multithreading
                    if (it == node2orbVec.end() || it->second.empty()) continue; //test debug test multithreading
                } // test debug test multithreading
                if (parindexVec_ref[n] < 0)
                    csize = sizecoeff;
                else
                    csize = sizecoeffW;

                // In the serial case we copy the coeff coeffBlock. In the mpi case coeffBlock is provided by the bank
                if (serial) {
                    int shift = sizecoeff - sizecoeffW; // to copy only wavelet part
                    if (parindexVec_ref[n] < 0) shift = 0;
                    DoubleMatrix coeffBlock(csize, it->second.size());
                    for (int j : it->second) { // loop over indices of the orbitals using this node
                        int orb_node_ix = orb2node[j][node_ix];
                        for (int k = 0; k < csize; k++) coeffBlock(k, orbVec.size()) = coeffVec[j][orb_node_ix][k + shift];
                        orbVec.push_back(j);
                    }
    // #pragma omp critical //Debug test debug (multithreading) seems to be ok now? 
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
            Fac[i] = BraKet[i].func_ptr->data.c1[l];
        }
        mrcpp::mpi::allreduce_vector(Fac, mrcpp::mpi::comm_wrk);
        for (int i = 0; i < N; i++) {
            if (!mrcpp::mpi::my_func(i)) continue;
            for (int j = 0; j < N; j++) { S(i, j) *= std::conj(Fac[i]) * Fac[j]; }
        }
        Stot += S;
    }
    return Stot;
}

/** @brief Compute the overlap matrix S_ij = <bra_i|ket_j>
 *
 *  Will take the conjugate of bra before integrating
 */
ComplexMatrix calc_overlap_matrix_cplx(CompFunctionVector &Bra, CompFunctionVector &Ket) {

    mrcpp::mpi::barrier(mrcpp::mpi::comm_wrk); // for consistent timings
    bool braisreal = false;
    bool ketisreal = false;
    for (int i = 0; i < Bra.size(); i++) if(Bra[i].isreal()) braisreal=true;
    for (int i = 0; i < Ket.size(); i++) if(Ket[i].isreal()) ketisreal=true;
    // copy real trees into complex ones, because it is much simpler than keeping tabs on 
    // which element of bra and ket is real 
    if (braisreal) {
        for (int i = 0; i < Bra.size(); i++) {
            if (!mrcpp::mpi::my_func(i)) continue;
            if (Bra[i].isreal()){
                for (int comp = 0; comp< Bra[i].Ncomp(); comp++) {
                    Bra[i].CompC[comp] = Bra[i].CompD[comp]->CopyTreeToComplex();
                    delete Bra[i].CompD[comp];
                    Bra[i].CompD[comp] = nullptr;
                }
                Bra[i].func_ptr->iscomplex = 1;
                Bra[i].func_ptr->isreal = 0; //just in case, to avoid defining the function as both real and complex, even though its real components are not empty
            }
        }
    }
    if (ketisreal) {
        for (int i = 0; i < Ket.size(); i++) {
            if (!mrcpp::mpi::my_func(i)) continue;
            if (Ket[i].isreal()){
                for (int comp = 0; comp< Ket[i].Ncomp(); comp++) {
                    Ket[i].CompC[comp] = Ket[i].CompD[comp]->CopyTreeToComplex();
                    delete Ket[i].CompD[comp];
                    Ket[i].CompD[comp] = nullptr;
                }
                Ket[i].func_ptr->iscomplex = 1;
                Ket[i].func_ptr->isreal = 0; //just in case, to avoid defining the function as both real and complex, even though its real components are not empty
            }
        }
    }

    std::shared_ptr<MultiResolutionAnalysis<3>> mra = Bra.vecMRA;

    int N = Bra.size();
    int M = Ket.size();
    ComplexMatrix Stot = ComplexMatrix::Zero(N, M);

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
    
    // Computing the number of components of the owned elements
    int Ncomponents = 1;
    for (int i = 0; i < N; i++) {
        if (!mrcpp::mpi::my_func(Bra[i])) continue;
        Ncomponents = std::max(Ncomponents, Bra[i].Ncomp()); //We assume Bra/Ket have the same number of components. It would anyway be problematic otherwise
    }
    for (int l = 0; l < Ncomponents; l++) {
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
        ComplexMatrix S = ComplexMatrix::Zero(N, M);
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
                Bra[j].complex(l).makeCoeffVector(coeffVecBra[j], indexVec, parindexVec, scalefac, max_ix, refTree);
                // make a map that gives j from indexVec
                int orb_node_ix = 0;
                for (int ix : indexVec) {
                    orb2nodeBra[j][ix] = orb_node_ix++;
                    if (ix < 0) continue;
                    node2orbVecBra[ix].push_back(j);
                }
            }
            for (int j = 0; j < M; j++) {
                // make vector with all coef pointers and their indices in the union grid
                Ket[j].complex(l).makeCoeffVector(coeffVecKet[j], indexVec, parindexVec, scalefac, max_ix, refTree);
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
            save_nodes(Bra, refTree, nodesBra, -1, l);
            save_nodes(Ket, refTree, nodesKet, -1, l);
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
                 
                int node_ix = indexVec_ref[n];      // SerialIx for this node in the reference tree
                std::map<int, std::vector<int>>::const_iterator itbra;
                std::map<int, std::vector<int>>::const_iterator itket;
                if (serial){ 
                    itbra = node2orbVecBra.find(node_ix); 
                    if (itbra == node2orbVecBra.end() || itbra->second.empty()) continue; 
                    itket = node2orbVecKet.find(node_ix); 
                    if (itket == node2orbVecKet.end() || itket->second.empty()) continue; 
                } 
                if (parindexVec_ref[n] < 0)
                    csize = sizecoeff;
                else
                    csize = sizecoeffW;
                if (serial) {
                    int shift = sizecoeff - sizecoeffW; // to copy only wavelet part
                    if (node2orbVecBra.count(node_ix) == 0) continue;
                    if (node2orbVecKet.count(node_ix) == 0) continue;
                    ComplexMatrix coeffBlockBra(csize, itbra->second.size()); 
                    ComplexMatrix coeffBlockKet(csize, itket->second.size()); 
                    if (parindexVec_ref[n] < 0) shift = 0;

                    for (int j : itbra->second) { 
                        int orb_node_ix = orb2nodeBra[j].at(node_ix); 
                        for (int k = 0; k < csize; k++) coeffBlockBra(k, orbVecBra.size()) = coeffVecBra[j][orb_node_ix][k + shift];
                        orbVecBra.push_back(j);
                    }
                    for (int j : itket->second) { 
                        int orb_node_ix = orb2nodeKet[j].at(node_ix); 
                        for (int k = 0; k < csize; k++) coeffBlockKet(k, orbVecKet.size()) = coeffVecKet[j][orb_node_ix][k + shift];
                        orbVecKet.push_back(j);
                    }
                    if (orbVecBra.size() > 0 and orbVecKet.size() > 0) {
                        ComplexMatrix S_temp(orbVecBra.size(), orbVecKet.size());
                        if (not conjMatBra[0] and not conjMatKet[0]) {
                            S_temp.noalias() = coeffBlockBra.transpose().conjugate() * coeffBlockKet;
                        } else if (conjMatBra[0] and not conjMatKet[0]) {
                            S_temp.noalias() = coeffBlockBra.transpose() * coeffBlockKet;
                        } else if (not conjMatBra[0] and conjMatKet[0]) {
                            S_temp.noalias() = coeffBlockBra.transpose() * coeffBlockKet.transpose();
                        } else if (conjMatBra[0] and conjMatKet[0]) {
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
                        if (not conjMatBra[0] and not conjMatKet[0]) {
                            S_temp.noalias() = coeffBlockBra.transpose().conjugate() * coeffBlockKet;
                        } else if (conjMatBra[0] and not conjMatKet[0]) {
                            S_temp.noalias() = coeffBlockBra.transpose() * coeffBlockKet;
                        } else if (not conjMatBra[0] and conjMatKet[0]) {
                            S_temp.noalias() = coeffBlockBra.transpose() * coeffBlockKet.transpose();
                        } else if (conjMatBra[0] and conjMatKet[0]) {
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
            FacBra[i] = Bra[i].func_ptr->data.c1[l];
        }
        for (int i = 0; i < M; i++) {
            if (!mrcpp::mpi::my_func(Ket[i])) continue;
            FacKet[i] = Ket[i].func_ptr->data.c1[l];
        }
        mrcpp::mpi::allreduce_vector(FacBra, mrcpp::mpi::comm_wrk);
        mrcpp::mpi::allreduce_vector(FacKet, mrcpp::mpi::comm_wrk);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) { S(i, j) *= std::conj(FacBra[i]) * FacKet[j]; }
        }
        Stot += S;
    }

    //THIS IS SEGMENTATION FAULT HELL! If for some reason one element, or worse, one component, 
    // of bra/ket is real but not the rest we simply nuke it for no reason
    // restore input
    // if (braisreal) {
    //     for (int i = 0; i < Bra.size(); i++) {
    //         for (int comp = 0; comp< Bra[k].Ncomp(); comp++) {
    //             delete Bra[i].CompC[comp];
    //             Bra[i].CompC[comp] = nullptr;
    //         }
    //         Bra[i].func_ptr->iscomplex = 0;
    //         Bra[i].func_ptr->isreal = 1;
    //     }
    // }
    // if (ketisreal) {
    //     for (int i = 0; i < Ket.size(); i++) {
    //         for (int comp = 0; comp< Ket[k].Ncomp(); comp++) {
    //             delete Ket[i].CompC[comp];
    //             Ket[i].CompC[comp] = nullptr;
    //         }
    //         Ket[i].func_ptr->iscomplex = 0;
    //         Ket[i].func_ptr->isreal = 1;
    //     }
    // }
    return Stot;
}

/** @brief Compute the overlap matrix S_ij = <bra_i|ket_j>
 *
 */
ComplexMatrix calc_overlap_matrix(CompFunctionVector &Bra, CompFunctionVector &Ket) { 
    int N = Bra.size();
    int M = Ket.size();

    bool bracomplex = false, ketcomplex = false;
    for (int i=0; i < N; i++) {
        if (!mrcpp::mpi::my_func(Bra[i])) continue;
        if (Bra[i].iscomplex()) bracomplex = true;
    }
    for (int i=0; i < M; i++) {
        if (!mrcpp::mpi::my_func(Ket[i])) continue;
        if (Ket[i].iscomplex()) ketcomplex = true;
    }
    bool complexcase = bracomplex or ketcomplex;
    complexcase = mrcpp::mpi::allreduce_max(complexcase ? 1 : 0, mrcpp::mpi::comm_wrk) > 0;
    if (complexcase) return calc_overlap_matrix_cplx(Bra, Ket);

    mrcpp::mpi::barrier(mrcpp::mpi::comm_wrk); // for consistent timings

    std::shared_ptr<MultiResolutionAnalysis<3>> mra = Bra.vecMRA;

    ComplexMatrix Stot = ComplexMatrix::Zero(N, M);

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

    // Computing the number of components of the owned elements
    int Ncomponents = 1;
    for (int i = 0; i < N; i++) {
        if (!mrcpp::mpi::my_func(Bra[i])) continue;
        Ncomponents = std::max(Ncomponents, Bra[i].Ncomp()); //we assume Bra and Ket to have the same number of component. It would be a problem if they didn't anyway.
    }
    // In the serial case we store the coeff pointers in coeffVec. In the mpi case the coeff are stored in the bank
    for (int l = 0; l < Ncomponents; l++) {
        // only used for serial case:
        std::vector<std::vector<double *>> coeffVecBra(N); //diabolus ex pointera 
        std::map<int, std::vector<int>> node2orbVecBra; // for each node index, gives a vector with the indices of the orbitals using this node
        std::vector<std::map<int, int>> orb2nodeBra(N); // for a given orbital and a given node, gives the node index in
                                                        // the orbital given the node index in the reference tree
        std::vector<std::vector<double *>> coeffVecKet(M);
        std::map<int, std::vector<int>> node2orbVecKet; // for each node index, gives a vector with the indices of the orbitals using this node
        std::vector<std::map<int, int>> orb2nodeKet(M); // for a given orbital and a given node, gives the node index in
                                                        // the orbital given the node index in the reference tree
        mrcpp::BankAccount nodesBra;
        mrcpp::BankAccount nodesKet;

        ComplexMatrix S = ComplexMatrix::Zero(N, M); // contribution to Stot from component k
        if (serial) {
            // 2) make list of all coefficients, and their reference indices
            // for different orbitals, indexVec will give the same index for the same node in space
            // TODO? : do not copy coefficients, but use directly the pointers
            // could OMP parallelize, but is fast anyway
            std::vector<int> parindexVec; // serialIx of the parent nodes
            std::vector<int> indexVec;    // serialIx of the nodes
            for (int j = 0; j < N; j++) {
                // make vector with all coef pointers and their indices in the union grid
                Bra[j].real(l).makeCoeffVector(coeffVecBra[j], indexVec, parindexVec, scalefac, max_ix, refTree);
                // make a map that gives j from indexVec
                int orb_node_ix = 0;
                for (int ix : indexVec) {
                    orb2nodeBra[j][ix] = orb_node_ix++;
                    if (ix < 0) continue;
                    node2orbVecBra[ix].push_back(j);
                }
            }
            for (int j = 0; j < M; j++) {
                // make vector with all coef pointers and their indices in the union grid
                Ket[j].real(l).makeCoeffVector(coeffVecKet[j], indexVec, parindexVec, scalefac, max_ix, refTree);
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
            save_nodes(Bra, refTree, nodesBra, -1, l);
            save_nodes(Ket, refTree, nodesKet, -1, l);
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
                // std::cout << "Processing node " << n + 1 << "/" << max_n << S(0,0) << std::endl;
                if (n % mrcpp::mpi::wrk_size != mrcpp::mpi::wrk_rank) continue;
                int csize;
                std::vector<int> orbVecBra; // identifies which Bra orbitals use this node
                std::vector<int> orbVecKet; // identifies which Ket orbitals use this node
                int node_ix = indexVec_ref[n];      // SerialIx for this node in the reference tree
                std::map<int, std::vector<int>>::const_iterator itbra; 
                std::map<int, std::vector<int>>::const_iterator itket; 
                if (serial){
                    itbra = node2orbVecBra.find(node_ix);
                    if (itbra == node2orbVecBra.end() || itbra->second.empty()) continue;
                    itket = node2orbVecKet.find(node_ix);
                    if (itket == node2orbVecKet.end() || itket->second.empty()) continue;
                } 
                if (parindexVec_ref[n] < 0)
                    csize = sizecoeff;
                else
                    csize = sizecoeffW;
                if (serial) {
                    int shift = sizecoeff - sizecoeffW; // to copy only wavelet part (no scaling function contribution)
                    if (node2orbVecBra.count(node_ix) == 0) continue;
                    if (node2orbVecKet.count(node_ix) == 0) continue;
                    DoubleMatrix coeffBlockBra(csize, itbra->second.size());
                    DoubleMatrix coeffBlockKet(csize, itket->second.size());
                    if (parindexVec_ref[n] < 0) shift = 0;
                    for (int j : itbra->second) { 
                        int orb_node_ix = orb2nodeBra[j].at(node_ix);
                        for (int k = 0; k < csize; k++) coeffBlockBra(k, orbVecBra.size()) = coeffVecBra[j][orb_node_ix][k + shift];
                        orbVecBra.push_back(j);
                    }
                    for (int j : itket->second) { 
                        int orb_node_ix = orb2nodeKet[j].at(node_ix);
                        for (int k = 0; k < csize; k++) coeffBlockKet(k, orbVecKet.size()) = coeffVecKet[j][orb_node_ix][k + shift];
                        orbVecKet.push_back(j);
                    }
    // #pragma omp critical //Debug test debug (multithreading) (copy fix from calc_overlap_matrix(braket), hope it works here as well)
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
                    for (int j = 0; j < M; j++) { 
                        S(i, j) += S_omp(i, j); 
                    }
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
            FacBra[i] = Bra[i].func_ptr->data.c1[l];
        }
        for (int i = 0; i < M; i++) {
            if (!mrcpp::mpi::my_func(Ket[i])) continue;
            FacKet[i] = Ket[i].func_ptr->data.c1[l];
        }
        mrcpp::mpi::allreduce_vector(FacBra, mrcpp::mpi::comm_wrk);
        mrcpp::mpi::allreduce_vector(FacKet, mrcpp::mpi::comm_wrk);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) { 
                S(i, j) *= std::conj(FacBra[i]) * FacKet[j]; 
            }
        }
        Stot += S; // accumulate contribution from this component
    }

    return Stot;
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

// /** @brief In-place element-wise addition of two CompFunctionVectors
//  * output[i] = out[i] + c[i]*inp[i]
//  *  @param out: receiving CompFunctionVector, will hold the result, but should not start out empty!
//  *  @param c: vector of coefficients 
//  *  @param inp: CompFunctionVector to be added
//  */
// void add(CompFunctionVector &out, ComplexVector c, CompFunctionVector &inp){
//     MSG_INFO("test install");
//     if (out.size()!=inp.size() or c.size()!=inp.size()) MSG_ABORT("Mismatched vector sizes for addition!");
//     for (int i=0; i < out.size(); i++){
//         //in-place addition
//         out[i].add(c[i], inp[i]);
//     }
// }

template void make_density(CompFunction<3> &out, CompFunction<3> &inp, double prec,  std::vector<bool> contrib = std::vector<bool>(4, true));
template ComplexDouble dot(const CompFunction<3> &bra, const CompFunction<3> &ket);
template void project(CompFunction<3> &out, RepresentableFunction<3, double> &f, double prec, int comp = 0);
template void project(CompFunction<3> &out, RepresentableFunction<3, ComplexDouble> &f, double prec, int comp = 0);
template void multiply(CompFunction<3> &out, CompFunction<3> inp_a, CompFunction<3> inp_b, double prec, bool absPrec, bool useMaxNorms, bool conjugate);
template void multiply(double prec, CompFunction<3> &out, double coef, CompFunction<3> inp_a, CompFunction<3> inp_b, int maxIter = -1, bool absPrec = false, bool useMaxNorms = false, bool conjugate = false);
template void multiply(CompFunction<3> &out, FunctionTree<3, double> &inp_a, RepresentableFunction<3, double> &f, double prec, int nrefine = 0, bool conjugate);
template void multiply(CompFunction<3> &out, FunctionTree<3, ComplexDouble> &inp_a, RepresentableFunction<3, ComplexDouble> &f, double prec, int nrefine = 0, bool conjugate);
template void multiply(CompFunction<3> &out, CompFunction<3> &inp_a, RepresentableFunction<3, double> &f, double prec, int nrefine = 0, bool conjugate);
template void multiply(CompFunction<3> &out, CompFunction<3> &inp_a, RepresentableFunction<3, ComplexDouble> &f, double prec, int nrefine = 0, bool conjugate);
template  void multiply(CompFunction<3> &out, CompFunction<3> inp_a, FunctionTree<3, double> &inp_b, double prec, bool absPrec = false, bool useMaxNorm = false, bool conjugate = false);
template  void multiply(CompFunction<3> &out, CompFunction<3> inp_a, FunctionTree<3, ComplexDouble> &inp_b, double prec, bool absPrec = false, bool useMaxNorm = false, bool conjugate = false);
template void CopyToComplex(CompFunction<3> &out, const CompFunction<3> &inp);
template void deep_copy(CompFunction<3> *out, const CompFunction<3> &inp);
template void deep_copy(CompFunction<3> &out, const CompFunction<3> &inp);
template void add(CompFunction<3> &out, ComplexDouble a, CompFunction<3> inp_a, ComplexDouble b, CompFunction<3> inp_b, double prec, bool conjugate);
template void linear_combination(CompFunction<3> &out, const std::vector<ComplexDouble> &c, std::vector<CompFunction<3>> &inp, double prec, bool conjugate);
template double node_norm_dot(CompFunction<3> bra, CompFunction<3> ket);
template void orthogonalize(double prec, CompFunction<3> &Bra, CompFunction<3> &Ket);
// template void make_density(CompFunction<3> &out, CompFunction<3> &inp, double prec);

} // namespace mrcpp
