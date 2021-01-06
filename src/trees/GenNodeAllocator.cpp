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

#include "GenNodeAllocator.h"
#include "GenNode.h"
#include "FunctionTree.h"
#include "utils/Printer.h"
#include "utils/mpi_utils.h"

namespace mrcpp {

/** SerialTree class constructor.
 * Allocate the root FunctionNodes and fill in the empty slots of rootBox.
 * Initializes rootNodes to represent the zero function and allocate their nodes.
 * NOTES:
 * Serial trees are made of projected nodes, and include gennodes and loose nodes separately.
 * All created (using class creator) Projected nodes or GenNodes are loose nodes.
 * Loose nodes have their coeff in serial Tree, but not the node part.
 * Projected nodes and GenNodes that are created by their creator, are detroyed by destructor ~ProjectedNode and
 * ~GenNode. Serial tree nodes are not using the destructors, but explicitely call to deallocNodes or deallocGenNodes
 * Gen nodes and loose nodes are not counted with MWTree->[in/de]crementNodeCount()
 */
template <int D>
GenNodeAllocator<D>::GenNodeAllocator(FunctionTree<D> *tree)
        : NodeAllocator<D>(tree, nullptr) {

    // Size for GenNodes chunks
    this->coeffsPerNode = this->tree_p->getKp1_d();       // One block
    println(10, "SizeNode Coeff (kB) " << this->coeffsPerNode * sizeof(double) / 1024);

    // 2 MB small for no waisting place, but large enough so that latency and overhead work is negligible
    int sizePerChunk = 2 * 1024 * 1024;
    if (D < 3) {
        // define rather from number of nodes per chunk
        this->maxNodesPerChunk = 64;
        sizePerChunk = this->maxNodesPerChunk * this->coeffsPerNode;
    } else {
        this->maxNodesPerChunk = (sizePerChunk / this->coeffsPerNode / sizeof(double) / 8) * 8;
    }

    // position just after last allocated node, i.e. where to put next node
    this->lastNode = this->sNodes;

    // make virtual table pointers
    auto *tmpGenNode = new GenNode<D>();
    this->cvptr_GenNode = *(char **)(tmpGenNode);
    delete tmpGenNode;

    MRCPP_INIT_OMP_LOCK();
}

template <int D> GenNodeAllocator<D>::~GenNodeAllocator() {
    for (int i = 0; i < this->nodeChunks.size(); i++) delete[](char *)(this->nodeChunks[i]);
    for (int i = 0; i < this->nodeCoeffChunks.size(); i++) delete[] this->nodeCoeffChunks[i];
    this->nodeStackStatus.clear();
    MRCPP_DESTROY_OMP_LOCK();
}

template <int D> void GenNodeAllocator<D>::allocRoots(MWTree<D> &tree) {
    NOT_REACHED_ABORT;
}

template <int D> void GenNodeAllocator<D>::allocChildrenNoCoeff(MWNode<D> &parent) {
    NOT_REACHED_ABORT;
}

template <int D> void GenNodeAllocator<D>::allocChildren(MWNode<D> &parent) {
    int sIx;
    double *coefs_p;
    // NB: serial tree MUST generate all children consecutively
    // all children must be generated at once if several threads are active
    int nChildren = parent.getTDim();
    GenNode<D> *child_p = this->allocNodes(nChildren, &sIx, &coefs_p);

    // position of first child
    parent.childSerialIx = sIx; // not used fro Gennodes?
    for (int cIdx = 0; cIdx < nChildren; cIdx++) {
        parent.children[cIdx] = child_p;

        *(char **)(child_p) = this->cvptr_GenNode;

        child_p->tree = parent.tree;
        child_p->parent = &parent;
        for (int i = 0; i < child_p->getTDim(); i++) { child_p->children[i] = nullptr; }

        child_p->maxSquareNorm = -1.0;
        child_p->maxWSquareNorm = -1.0;

        child_p->nodeIndex = parent.getNodeIndex().child(cIdx);
        child_p->hilbertPath = HilbertPath<D>(parent.getHilbertPath(), cIdx);

        child_p->n_coefs = this->coeffsPerNode;
        child_p->coefs = coefs_p;

        child_p->lockX = 0;
        child_p->serialIx = sIx;
        child_p->parentSerialIx = parent.serialIx;
        child_p->childSerialIx = -1;

        child_p->status = 0;

        child_p->clearNorms();
        child_p->setIsLeafNode();
        child_p->setIsAllocated();
        child_p->clearHasCoefs();
        child_p->setIsGenNode();

        child_p->tree->incrementGenNodeCount();

        sIx++;
        child_p++;
        coefs_p += this->coeffsPerNode;
    }
}

// return pointer to the last active node or NULL if failed
template <int D> GenNode<D> *GenNodeAllocator<D>::allocNodes(int nAlloc, int *serialIx, double **coefs_p) {
    MRCPP_SET_OMP_LOCK();
    *serialIx = this->nNodes;
    int chunkIx = *serialIx % (this->maxNodesPerChunk);

    // Not necessarily wrong, but new:
    assert(nAlloc == (1 << D));

    if (chunkIx == 0 or chunkIx + nAlloc > this->maxNodesPerChunk) {
        // start on new chunk
        // we want nodes allocated simultaneously to be allocated on the same chunk.
        // possibly jump over the last nodes from the old chunk
        this->nNodes = this->maxNodesPerChunk * ((this->nNodes + nAlloc - 1) / this->maxNodesPerChunk);

        int chunk = this->nNodes / this->maxNodesPerChunk; // find the right chunk

        // careful: nodeChunks.size() is an unsigned int
        if (chunk + 1 > this->nodeChunks.size()) {
            // need to allocate new chunk
            this->sNodes = (GenNode<D> *)new char[this->maxNodesPerChunk * sizeof(GenNode<D>)];
            for (int i = 0; i < this->maxNodesPerChunk; i++) {
                this->sNodes[i].serialIx = -1;
                this->sNodes[i].parentSerialIx = -1;
                this->sNodes[i].childSerialIx = -1;
            }
            this->nodeChunks.push_back(this->sNodes);
            auto *sNodesCoeff = new double[this->coeffsPerNode * this->maxNodesPerChunk];
            this->nodeCoeffChunks.push_back(sNodesCoeff);
            // allocate new chunk in nodeStackStatus
            int oldsize = this->nodeStackStatus.size();
            int newsize = oldsize + this->maxNodesPerChunk;
            for (int i = oldsize; i < newsize; i++) this->nodeStackStatus.push_back(0);
            this->maxNodes = newsize;

            if (chunk % 100 == 99 and D == 3)
                println(10,
                        "\n number of GenNodes " << this->nNodes << ",number of GenNodechunks now "
                                                 << this->nodeChunks.size() << ", total size coeff  (MB) "
                                                 << (this->nNodes / 1024) * this->coeffsPerNode / 128);
        }
        this->lastNode = this->nodeChunks[chunk] + this->nNodes % (this->maxNodesPerChunk);
        *serialIx = this->nNodes;
        chunkIx = *serialIx % (this->maxNodesPerChunk);
    }
    assert((this->nNodes + nAlloc - 1) / this->maxNodesPerChunk < this->nodeChunks.size());

    GenNode<D> *newNode = this->lastNode;
    GenNode<D> *newNode_cp = newNode;

    int chunk = this->nNodes / this->maxNodesPerChunk; // find the right chunk
    *coefs_p = this->nodeCoeffChunks[chunk] + chunkIx * this->coeffsPerNode;

    for (int i = 0; i < nAlloc; i++) {
        newNode_cp->serialIx = *serialIx + i; // Until overwritten!
        if (this->nodeStackStatus[*serialIx + i] != 0)
            println(0, *serialIx + i << " NodeStackStatus: not available " << this->nodeStackStatus[*serialIx + i]);
        this->nodeStackStatus[*serialIx + i] = 1;
        newNode_cp++;
    }
    this->nNodes += nAlloc;
    this->lastNode += nAlloc;

    MRCPP_UNSET_OMP_LOCK();
    return newNode;
}

template <int D> void GenNodeAllocator<D>::deallocNodes(int serialIx) {
    MRCPP_SET_OMP_LOCK();
    if (this->nNodes < 0) {
        println(0, "minNodes exceeded " << this->nNodes);
        this->nNodes++;
    }
    this->nodeStackStatus[serialIx] = 0; // mark as available
    if (serialIx == this->nNodes - 1) {  // top of stack
        int topStack = this->nNodes;
        while (this->nodeStackStatus[topStack - 1] == 0) {
            topStack--;
            if (topStack < 1) {
                // remove all the NodeChunks once there are noe more GenNodes
                this->deallocNodeChunks();
                break;
            }
        }
        this->nNodes = topStack; // move top of stack
        // has to redefine lastGenNode
        int chunk = this->nNodes / this->maxNodesPerChunk; // find the right chunk
        this->lastNode = this->nodeChunks[chunk] + this->nNodes % (this->maxNodesPerChunk);
    }
    MRCPP_UNSET_OMP_LOCK();
}

template <int D> void GenNodeAllocator<D>::deallocNodeChunks() {
    for (int i = 0; i < this->nodeCoeffChunks.size(); i++) delete[] this->nodeCoeffChunks[i];
    for (int i = 0; i < this->nodeChunks.size(); i++) delete[](char *)(this->nodeChunks[i]);
    this->nodeCoeffChunks.clear();
    this->nodeChunks.clear();
    this->nodeStackStatus.clear();
}

template class GenNodeAllocator<1>;
template class GenNodeAllocator<2>;
template class GenNodeAllocator<3>;
} // namespace mrcpp
