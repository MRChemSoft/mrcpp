/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2020 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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

#include <fstream>

#include "FunctionNode.h"
#include "FunctionTree.h"
#include "HilbertIterator.h"
#include "ProjectedNode.h"
#include "SerialFunctionTree.h"
#include "utils/Printer.h"
#include "utils/Timer.h"
#include "utils/mpi_utils.h"
#include "utils/periodic_utils.h"

using namespace Eigen;

namespace mrcpp {

/** @returns New FunctionTree object
 *
 *  @param[in] mra: Which MRA the function is defined
 *  @param[in] sh_mem: Pointer to MPI shared memory block
 *
 *  @details Constructs an uninitialized tree, containing only empty root nodes.
 *  If a shared memory pointer is provided the tree will be allocated in this
 *  shared memory window, otherwise it will be local to each MPI process.
 */
template <int D>
FunctionTree<D>::FunctionTree(const MultiResolutionAnalysis<D> &mra, SharedMemory *sh_mem)
        : MWTree<D>(mra)
        , RepresentableFunction<D>(mra.getWorldBox().getLowerBounds().data(),
                                   mra.getWorldBox().getUpperBounds().data()) {
    this->serialTree_p = new SerialFunctionTree<D>(this, sh_mem);
    this->serialTree_p->allocRoots(*this);
    this->resetEndNodeTable();
}

// FunctionTree destructor
template <int D> FunctionTree<D>::~FunctionTree() {
    for (int i = 0; i < this->rootBox.size(); i++) {
        MWNode<D> &root = this->getRootMWNode(i);
        root.deleteChildren();
        root.dealloc();
        this->rootBox.clearNode(i);
    }
    delete this->serialTree_p;
}

/** @brief Remove all nodes in the tree
 *
 * @details Leaves the tree inn the same state as after construction, i.e.
 * undefined function containing only root nodes without coefficients.
 * The assigned memory (nodeChunks in SerialTree) is NOT released,
 * but is immediately available to the new function.
 */
template <int D> void FunctionTree<D>::clear() {
    for (int i = 0; i < this->rootBox.size(); i++) {
        MWNode<D> &root = this->getRootMWNode(i);
        root.deleteChildren();
        root.clearHasCoefs();
        root.clearNorms();
    }
    this->resetEndNodeTable();
    this->clearSquareNorm();
    this->getSerialFunctionTree()->clear(this->rootBox.size());
}

/** @brief Write the tree structure to disk, for later use
 * @param[in] file: File name, will get ".tree" extension
 */
template <int D> void FunctionTree<D>::saveTree(const std::string &file) {
    // This is basically a copy of MPI send_tree
    Timer t1;
    std::stringstream fname;
    fname << file << ".tree";

    std::fstream f;
    f.open(fname.str(), std::ios::out | std::ios::binary);
    if (not f.is_open()) MSG_ERROR("Unable to open file");

    this->deleteGenerated();
    SerialFunctionTree<D> &sTree = *this->getSerialFunctionTree();

    // Write size of tree
    int nChunks = sTree.getNChunksUsed();
    f.write((char *)&nChunks, sizeof(int));

    // Write tree data, chunk by chunk
    int count = 1;
    for (int iChunk = 0; iChunk < nChunks; iChunk++) {
        count = sTree.maxNodesPerChunk * sizeof(ProjectedNode<D>);
        f.write((char *)sTree.nodeChunks[iChunk], count);
        count = sTree.sizeNodeCoeff * sTree.maxNodesPerChunk;
        f.write((char *)sTree.nodeCoeffChunks[iChunk], count * sizeof(double));
    }
    f.close();
    print::time(10, "Time write", t1);
}

/** @brief Read a previously stored tree structure from disk
 * @param[in] file: File name, will get ".tree" extension
 * @note This tree must have the exact same MRA the one that was saved
 */
template <int D> void FunctionTree<D>::loadTree(const std::string &file) {
    // This is basically a copy of MPI recv_tree
    Timer t1;
    std::stringstream fname;
    fname << file << ".tree";

    std::fstream f;
    f.open(fname.str(), std::ios::in | std::ios::binary);
    if (not f.is_open()) MSG_ERROR("Unable to open file");

    // Read size of tree
    int nChunks;
    f.read((char *)&nChunks, sizeof(int));
    SerialFunctionTree<D> &sTree = *this->getSerialFunctionTree();

    // Read tree data, chunk by chunk
    int count = 1;
    for (int iChunk = 0; iChunk < nChunks; iChunk++) {
        if (iChunk < sTree.nodeChunks.size()) {
            sTree.sNodes = sTree.nodeChunks[iChunk];
        } else {
            double *sNodesCoeff;
            if (sTree.isShared()) {
                // for coefficients, take from the shared memory block
                SharedMemory *shMem = sTree.getMemory();
                sNodesCoeff = shMem->sh_end_ptr;
                shMem->sh_end_ptr += (sTree.sizeNodeCoeff * sTree.maxNodesPerChunk);
                // may increase size dynamically in the future
                if (shMem->sh_max_ptr < shMem->sh_end_ptr) MSG_ABORT("Shared block too small");
            } else {
                sNodesCoeff = new double[sTree.sizeNodeCoeff * sTree.maxNodesPerChunk];
            }
            sTree.nodeCoeffChunks.push_back(sNodesCoeff);
            sTree.sNodes = (ProjectedNode<D> *)new char[sTree.maxNodesPerChunk * sizeof(ProjectedNode<D>)];
            sTree.nodeChunks.push_back(sTree.sNodes);
        }
        count = sTree.maxNodesPerChunk * sizeof(ProjectedNode<D>);
        f.read((char *)sTree.nodeChunks[iChunk], count);
        count = sTree.sizeNodeCoeff * sTree.maxNodesPerChunk;
        f.read((char *)sTree.nodeCoeffChunks[iChunk], count * sizeof(double));
    }
    f.close();
    print::time(10, "Time read tree", t1);

    Timer t2;
    sTree.rewritePointers();
    print::time(10, "Time rewrite pointers", t2);
}

/** @returns Integral of the function over the entire computational domain */
template <int D> double FunctionTree<D>::integrate() const {

    double result = 0.0;
    for (int i = 0; i < this->rootBox.size(); i++) {
        const FunctionNode<D> &fNode = getRootFuncNode(i);
        result += fNode.integrate();
    }

    // Handle potential scaling
    auto sf = this->getMRA().getWorldBox().getScalingFactor();
    auto jacobian = 1.0;
    for (const auto &sf_i : sf) jacobian *= std::sqrt(sf_i);
    // Square root of scaling factor in each diection. The seemingly missing
    // multiplication by the square root of sf_i is included in the basis

    return jacobian * result;
}

/** @returns Function value in a point, out of bounds returns zero
 *
 * @param[in] r: Cartesian coordinate
 *
 * @note This will only evaluate the _scaling_ part of the
 *       leaf nodes in the tree, which means that the function
 *       values will not be fully accurate. If you want to include
 *       also the _final_ wavelet part you'll have to manually extend
 *       the MW grid by one level before evaluating, using
 *       `mrcpp::refine_grid(tree, 1)`
 *       This is done to allow a fast and const function evaluation
 *       that can be done in OMP parallel.
 */
template <int D> double FunctionTree<D>::evalf(const Coord<D> &r) const {
    // Handle potential scaling
    const auto sf = this->getMRA().getWorldBox().getScalingFactor();
    auto arg = r;
    for (auto i = 0; i < D; i++) arg[i] = arg[i] / sf[i];

    // Adjust for scaling factor included in basis
    auto coef = 1.0;
    for (const auto &fac : sf) coef /= std::sqrt(fac);

    // The 1.0 appearing in the if tests comes from the period is
    // always 1.0 from the point of view of this function.
    if (this->getRootBox().isPeriodic()) { periodic::coord_manipulation<D>(arg, this->getRootBox().getPeriodic()); }

    // Function is zero outside the domain
    if (this->outOfBounds(arg)) return 0.0;

    const MWNode<D> &mw_node = this->getNodeOrEndNode(arg);
    auto &f_node = static_cast<const FunctionNode<D> &>(mw_node);
    auto result = f_node.evalScaling(arg);
    return coef * result;
}

/** @brief In-place square of MW function representations, fixed grid
 *
 * @details The leaf node point values of the function will be in-place
 * squared, no grid refinement.
 *
 */
template <int D> void FunctionTree<D>::square() {
    if (this->getNGenNodes() != 0) MSG_ABORT("GenNodes not cleared");

#pragma omp parallel num_threads(mrcpp_get_num_threads())
    {
        int nNodes = this->getNEndNodes();
        int nCoefs = this->getTDim() * this->getKp1_d();
#pragma omp for schedule(guided)
        for (int n = 0; n < nNodes; n++) {
            MWNode<D> &node = *this->endNodeTable[n];
            node.mwTransform(Reconstruction);
            node.cvTransform(Forward);
            double *coefs = node.getCoefs();
            for (int i = 0; i < nCoefs; i++) { coefs[i] *= coefs[i]; }
            node.cvTransform(Backward);
            node.mwTransform(Compression);
            node.calcNorms();
        }
    }
    this->mwTransform(BottomUp);
    this->calcSquareNorm();
}

/** @brief In-place power of MW function representations, fixed grid
 *
 * @param[in] p: Numerical power
 *
 * @details The leaf node point values of the function will be in-place raised
 * to the given power, no grid refinement.
 *
 */
template <int D> void FunctionTree<D>::power(double p) {
    if (this->getNGenNodes() != 0) MSG_ABORT("GenNodes not cleared");

#pragma omp parallel num_threads(mrcpp_get_num_threads())
    {
        int nNodes = this->getNEndNodes();
        int nCoefs = this->getTDim() * this->getKp1_d();
#pragma omp for schedule(guided)
        for (int n = 0; n < nNodes; n++) {
            MWNode<D> &node = *this->endNodeTable[n];
            node.mwTransform(Reconstruction);
            node.cvTransform(Forward);
            double *coefs = node.getCoefs();
            for (int i = 0; i < nCoefs; i++) { coefs[i] = std::pow(coefs[i], p); }
            node.cvTransform(Backward);
            node.mwTransform(Compression);
            node.calcNorms();
        }
    }
    this->mwTransform(BottomUp);
    this->calcSquareNorm();
}

/** @brief In-place multiplication by a scalar, fixed grid
 *
 * @param[in] c: Scalar coefficient
 *
 * @details The leaf node point values of the function will be
 * in-place multiplied by the given coefficient, no grid refinement.
 *
 */
template <int D> void FunctionTree<D>::rescale(double c) {
    if (this->getNGenNodes() != 0) MSG_ABORT("GenNodes not cleared");
#pragma omp parallel firstprivate(c) num_threads(mrcpp_get_num_threads())
    {
        int nNodes = this->getNEndNodes();
        int nCoefs = this->getTDim() * this->getKp1_d();
#pragma omp for schedule(guided)
        for (int i = 0; i < nNodes; i++) {
            MWNode<D> &node = *this->endNodeTable[i];
            if (not node.hasCoefs()) MSG_ABORT("No coefs");
            double *coefs = node.getCoefs();
            for (int j = 0; j < nCoefs; j++) { coefs[j] *= c; }
            node.calcNorms();
        }
    }
    this->mwTransform(BottomUp);
    this->calcSquareNorm();
}

/** @brief In-place rescaling by a function norm \f$ ||f||^{-1} \f$, fixed grid */
template <int D> void FunctionTree<D>::normalize() {
    if (this->getNGenNodes() != 0) MSG_ABORT("GenNodes not cleared");
    double sq_norm = this->getSquareNorm();
    if (sq_norm < 0.0) MSG_ERROR("Normalizing uninitialized function");
    this->rescale(1.0 / std::sqrt(sq_norm));
}

/** @brief In-place addition with MW function representations, fixed grid
 *
 * @param[in] c: Numerical coefficient of input function
 * @param[in] inp: Input function to add
 *
 * @details The input function will be added in-place on the current grid of
 * the function, i.e. no further grid refinement.
 *
 */
template <int D> void FunctionTree<D>::add(double c, FunctionTree<D> &inp) {
    if (this->getMRA() != inp.getMRA()) MSG_ABORT("Incompatible MRA");
    if (this->getNGenNodes() != 0) MSG_ABORT("GenNodes not cleared");
#pragma omp parallel firstprivate(c) shared(inp) num_threads(mrcpp_get_num_threads())
    {
        int nNodes = this->getNEndNodes();
#pragma omp for schedule(guided)
        for (int n = 0; n < nNodes; n++) {
            MWNode<D> &out_node = *this->endNodeTable[n];
            MWNode<D> &inp_node = inp.getNode(out_node.getNodeIndex());
            double *out_coefs = out_node.getCoefs();
            const double *inp_coefs = inp_node.getCoefs();
            for (int i = 0; i < inp_node.getNCoefs(); i++) { out_coefs[i] += c * inp_coefs[i]; }
            out_node.calcNorms();
        }
    }
    this->mwTransform(BottomUp);
    this->calcSquareNorm();
    inp.deleteGenerated();
}
/** @brief In-place addition of absolute values of MW function representations
 *
 * @param[in] c Numerical coefficient of input function
 * @param[in] inp Input function to add
 *
 * The absolute value of input function will be added in-place on the current grid of the output
 * function, i.e. no further grid refinement.
 *
 */
template <int D> void FunctionTree<D>::absadd(double c, FunctionTree<D> &inp) {
    if (this->getNGenNodes() != 0) MSG_ABORT("GenNodes not cleared");
#pragma omp parallel firstprivate(c) shared(inp) num_threads(mrcpp_get_num_threads())
    {
        int nNodes = this->getNEndNodes();
#pragma omp for schedule(guided)
        for (int n = 0; n < nNodes; n++) {
            MWNode<D> &out_node = *this->endNodeTable[n];
            MWNode<D> inp_node = inp.getNode(out_node.getNodeIndex()); // Full copy
            out_node.mwTransform(Reconstruction);
            out_node.cvTransform(Forward);
            inp_node.mwTransform(Reconstruction);
            inp_node.cvTransform(Forward);
            double *out_coefs = out_node.getCoefs();
            const double *inp_coefs = inp_node.getCoefs();
            for (int i = 0; i < inp_node.getNCoefs(); i++) { out_coefs[i] = abs(out_coefs[i]) + c * abs(inp_coefs[i]); }
            out_node.cvTransform(Backward);
            out_node.mwTransform(Compression);
            out_node.calcNorms();
        }
    }
    this->mwTransform(BottomUp);
    this->calcSquareNorm();
    inp.deleteGenerated();
}

/** @brief In-place multiplication with MW function representations, fixed grid
 *
 * @param[in] c: Numerical coefficient of input function
 * @param[in] inp: Input function to multiply
 *
 * @details The input function will be multiplied in-place on the current grid
 * of the function, i.e. no further grid refinement.
 *
 */
template <int D> void FunctionTree<D>::multiply(double c, FunctionTree<D> &inp) {
    if (this->getMRA() != inp.getMRA()) MSG_ABORT("Incompatible MRA");
    if (this->getNGenNodes() != 0) MSG_ABORT("GenNodes not cleared");
#pragma omp parallel firstprivate(c) shared(inp) num_threads(mrcpp_get_num_threads())
    {
        int nNodes = this->getNEndNodes();
#pragma omp for schedule(guided)
        for (int n = 0; n < nNodes; n++) {
            MWNode<D> &out_node = *this->endNodeTable[n];
            MWNode<D> inp_node = inp.getNode(out_node.getNodeIndex()); // Full copy
            out_node.mwTransform(Reconstruction);
            out_node.cvTransform(Forward);
            inp_node.mwTransform(Reconstruction);
            inp_node.cvTransform(Forward);
            double *out_coefs = out_node.getCoefs();
            const double *inp_coefs = inp_node.getCoefs();
            for (int i = 0; i < inp_node.getNCoefs(); i++) { out_coefs[i] *= c * inp_coefs[i]; }
            out_node.cvTransform(Backward);
            out_node.mwTransform(Compression);
            out_node.calcNorms();
        }
    }
    this->mwTransform(BottomUp);
    this->calcSquareNorm();
    inp.deleteGenerated();
}

/** @brief In-place mapping with a predefined function f(x), fixed grid
 *
 * @param[in] fmap: mapping function
 *
 * @details The input function will be mapped in-place on the current grid
 * of the function, i.e. no further grid refinement.
 *
 */
template <int D> void FunctionTree<D>::map(FMap fmap) {
    if (this->getNGenNodes() != 0) MSG_ABORT("GenNodes not cleared");
    {
        int nNodes = this->getNEndNodes();
#pragma omp parallel for schedule(guided) num_threads(mrcpp_get_num_threads())
        for (int n = 0; n < nNodes; n++) {
            MWNode<D> &node = *this->endNodeTable[n];
            node.mwTransform(Reconstruction);
            node.cvTransform(Forward);
            double *coefs = node.getCoefs();
            for (int i = 0; i < node.getNCoefs(); i++) { coefs[i] = fmap(coefs[i]); }
            node.cvTransform(Backward);
            node.mwTransform(Compression);
            node.calcNorms();
        }
    }
    this->mwTransform(BottomUp);
    this->calcSquareNorm();
}

template <int D> int FunctionTree<D>::getNChunks() {
    return this->getSerialFunctionTree()->getNChunks();
}

template <int D> int FunctionTree<D>::getNChunksUsed() {
    return this->getSerialFunctionTree()->getNChunksUsed();
}

template <int D> void FunctionTree<D>::getEndValues(VectorXd &data) {
    if (this->getNGenNodes() != 0) MSG_ABORT("GenNodes not cleared");
    int nNodes = this->getNEndNodes();
    int nCoefs = this->getTDim() * this->getKp1_d();
    data = VectorXd::Zero(nNodes * nCoefs);
    for (int n = 0; n < nNodes; n++) {
        MWNode<D> &node = getEndFuncNode(n);
        node.mwTransform(Reconstruction);
        node.cvTransform(Forward);
        const double *c = node.getCoefs();
        for (int i = 0; i < nCoefs; i++) { data(n * nCoefs + i) = c[i]; }
        node.cvTransform(Backward);
        node.mwTransform(Compression);
    }
}

template <int D> void FunctionTree<D>::setEndValues(VectorXd &data) {
    if (this->getNGenNodes() != 0) MSG_ABORT("GenNodes not cleared");
    int nNodes = this->getNEndNodes();
    int nCoefs = this->getTDim() * this->getKp1_d();
    for (int i = 0; i < nNodes; i++) {
        MWNode<D> &node = getEndFuncNode(i);
        const double *c = data.segment(i * nCoefs, nCoefs).data();
        node.setCoefBlock(0, nCoefs, c);
        node.cvTransform(Backward);
        node.mwTransform(Compression);
        node.setHasCoefs();
        node.calcNorms();
    }
    this->mwTransform(BottomUp);
    this->calcSquareNorm();
}

template <int D> std::ostream &FunctionTree<D>::print(std::ostream &o) {
    o << std::endl << "*FunctionTree: " << this->name << std::endl;
    return MWTree<D>::print(o);
}

template <int D> void FunctionTree<D>::printSerialIndices() {
    SerialFunctionTree<D> &sTree = *this->getSerialFunctionTree();
    int n = 0;
    for (int iChunk = 0; iChunk < sTree.getNChunks(); iChunk++) {
        int iShift = iChunk * sTree.maxNodesPerChunk;
        printout(0, "new chunk \n");
        for (int i = 0; i < sTree.maxNodesPerChunk; i++) {
            int status = sTree.nodeStackStatus[iShift + i];
            int sIdx = sTree.nodeChunks[iChunk][i].serialIx;
            int pIdx = sTree.nodeChunks[iChunk][i].parentSerialIx;
            int cIdx = sTree.nodeChunks[iChunk][i].childSerialIx;
            printout(0, std::setw(4) << n++);
            printout(0, std::setw(4) << status);
            printout(0, std::setw(6) << sIdx);
            printout(0, std::setw(6) << pIdx);
            printout(0, std::setw(6) << cIdx << "   ");
            if (status == 1) printout(0, sTree.nodeChunks[iChunk][i].getNodeIndex());
            printout(0, "\n");
        }
    }
}

/** @brief Reduce the precision of the tree by deleting nodes
 *
 * @param prec: New precision criterion
 * @param splitFac: Splitting factor: 1, 2 or 3
 * @param absPrec: Use absolute precision
 *
 * @details This will run the tree building algorithm in "reverse", starting
 * from the leaf nodes, and perform split checks on each node based on the given
 * precision and the local wavelet norm.
 *
 * @note The splitting factor appears in the threshold for the wavelet norm as
 * \f$ ||w|| < 2^{-sn/2} ||f|| \epsilon \f$. In principal, `s` should be equal
 * to the dimension; in practice, it is set to `s=1`.
 */
template <int D> int FunctionTree<D>::crop(double prec, double splitFac, bool absPrec) {

    for (int i = 0; i < this->rootBox.size(); i++) {
        MWNode<D> &root = this->getRootMWNode(i);
        root.crop(prec, splitFac, absPrec);
    }
    int nChunks = this->getSerialFunctionTree()->shrinkChunks();
    this->resetEndNodeTable();
    this->calcSquareNorm();
    return nChunks;
}

template class FunctionTree<1>;
template class FunctionTree<2>;
template class FunctionTree<3>;

} // namespace mrcpp
