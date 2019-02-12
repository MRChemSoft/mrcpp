/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2019 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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

/**
 *  \date 2016
 *          CTCC, University of Tromsø
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

using namespace Eigen;

namespace mrcpp {

/** FunctionTree constructor for Serial Tree.
 * */
template <int D>
FunctionTree<D>::FunctionTree(const MultiResolutionAnalysis<D> &mra, SharedMemory *sh_mem)
        : MWTree<D>(mra) {
    this->serialTree_p = new SerialFunctionTree<D>(this, sh_mem);
    this->serialTree_p->allocRoots(*this);
    this->resetEndNodeTable();
}

//** FunctionTree destructor. */
template <int D> FunctionTree<D>::~FunctionTree() {
    for (int i = 0; i < this->rootBox.size(); i++) {
        MWNode<D> &root = this->getRootMWNode(i);
        root.deleteChildren();
        root.dealloc();
        this->rootBox.clearNode(i);
    }
    delete this->serialTree_p;
}

/** Leaves the tree inn the same state as after construction, e.i.
 * undefined function containing only root nodes without coefficients.
 * The assigned memory (nodeChunks in SerialTree) is NOT released,
 * but is immediately available to the new function. */
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

/** Write the tree structure to disk, for later use.
 * Argument file name will get a ".tree" file extension. */
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

    t1.stop();
    Printer::printTime(10, "Time write", t1);
}

/** Read a previously stored tree structure from disk.
 * Argument file name will get a ".tree" file extension. */
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
                if (shMem->sh_max_ptr < shMem->sh_end_ptr) { MSG_FATAL("Shared block too small"); }
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

    t1.stop();
    Printer::printTime(10, "Time read tree", t1);

    Timer t2;
    sTree.rewritePointers();
    t2.stop();
    Printer::printTime(10, "Time rewrite pointers", t2);
}

template <int D> double FunctionTree<D>::integrate() const {

    double result = 0.0;
    for (int i = 0; i < this->rootBox.size(); i++) {
        const FunctionNode<D> &fNode = getRootFuncNode(i);
        result += fNode.integrate();
    }

    // Hande potential scaling
    auto sf = this->getMRA().getWorldBox().getScalingFactor();
    auto jacobian = 1.0;
    for (const auto &sf_i : sf) jacobian *= std::sqrt(sf_i);
    // Square root of scaling factor in each diection. The seemingly missing
    // multiplication by the square root of sf_i is included in the basis

    return jacobian * result;
}

template <int D> double FunctionTree<D>::evalf(const Coord<D> &r) {

    // Hande potential scaling
    const auto sf = this->getMRA().getWorldBox().getScalingFactor();
    auto arg = r;
    for (auto i = 0; i < D; i++) arg[i] = arg[i] / sf[i];

    // Adjust for scaling factor included in basis
    auto coef = 1.0;
    for (const auto &fac : sf) coef /= std::sqrt(fac);

    MWNode<D> &mr_node = this->getNodeOrEndNode(arg);
    auto &f_node = static_cast<FunctionNode<D> &>(mr_node);
    auto result = f_node.evalf(arg);
    this->deleteGenerated();
    return coef * result;
}
/** @brief In-place square of function
 *
 * The leaf node point values of the output function will be in-place
 * squared, no grid refinement.
 *
 */
template <int D> void FunctionTree<D>::square() {
    if (this->getNGenNodes() != 0) MSG_FATAL("GenNodes not cleared");

#pragma omp parallel
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

/** @brief In-place raise to given power
 *
 * @param[in] c Numerical power
 *
 * The leaf node point values of the output function will be in-place
 * raised to the given power, no grid refinement.
 *
 */
template <int D> void FunctionTree<D>::power(double p) {
    if (this->getNGenNodes() != 0) MSG_FATAL("GenNodes not cleared");

#pragma omp parallel
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

/** @brief In-place multiplication by a scalar
 *
 * @param[in] c Scalar coefficient
 *
 * The leaf node point values of the output function will be in-place
 * multiplied by the given coefficient, no grid refinement.
 *
 */
template <int D> void FunctionTree<D>::rescale(double c) {
    if (this->getNGenNodes() != 0) MSG_FATAL("GenNodes not cleared");
#pragma omp parallel firstprivate(c)
    {
        int nNodes = this->getNEndNodes();
        int nCoefs = this->getTDim() * this->getKp1_d();
#pragma omp for schedule(guided)
        for (int i = 0; i < nNodes; i++) {
            MWNode<D> &node = *this->endNodeTable[i];
            if (not node.hasCoefs()) MSG_FATAL("No coefs");
            double *coefs = node.getCoefs();
            for (int j = 0; j < nCoefs; j++) { coefs[j] *= c; }
            node.calcNorms();
        }
    }
    this->mwTransform(BottomUp);
    this->calcSquareNorm();
}

template <int D> void FunctionTree<D>::normalize() {
    if (this->getNGenNodes() != 0) MSG_FATAL("GenNodes not cleared");
    double sq_norm = this->getSquareNorm();
    if (sq_norm < 0.0) MSG_ERROR("Normalizing uninitialized function");
    this->rescale(1.0 / std::sqrt(sq_norm));
}

/** @brief In-place addition of MW function representations
 *
 * @param[in] c Numerical coefficient of input function
 * @param[in] inp Input function to add
 *
 * The input function will be added in-place on the current grid of the output
 * function, i.e. no further grid refinement.
 *
 */
template <int D> void FunctionTree<D>::add(double c, FunctionTree<D> &inp) {
    if (this->getNGenNodes() != 0) MSG_FATAL("GenNodes not cleared");
#pragma omp parallel firstprivate(c), shared(inp)
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

/** @brief In-place multiplication of MW function representations
 *
 * @param[in] c Numerical coefficient of input function
 * @param[in] inp Input function to multiply
 *
 * The input function will be multiplied in-place on the current grid of the
 * output function, i.e. no further grid refinement.
 *
 */
template <int D> void FunctionTree<D>::multiply(double c, FunctionTree<D> &inp) {
    if (this->getNGenNodes() != 0) MSG_FATAL("GenNodes not cleared");
#pragma omp parallel firstprivate(c), shared(inp)
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

template <int D> int FunctionTree<D>::getNChunks() {
    return this->getSerialFunctionTree()->getNChunks();
}

template <int D> int FunctionTree<D>::getNChunksUsed() {
    return this->getSerialFunctionTree()->getNChunksUsed();
}

template <int D> void FunctionTree<D>::getEndValues(VectorXd &data) {
    if (this->getNGenNodes() != 0) MSG_FATAL("GenNodes not cleared");
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
    if (this->getNGenNodes() != 0) MSG_FATAL("GenNodes not cleared");
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

/** Reduce the accuracy of the tree by deleting nodes
 * which have a higher precision than the requested precison.
 * By default, the relative precision of the tree is used. */
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
