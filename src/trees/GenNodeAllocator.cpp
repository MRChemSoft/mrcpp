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
#include "FunctionTree.h"
#include "FunctionNode.h"
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
 * Projected nodes and GenNodes that are created by their creator, are detroyed by destructor ~FunctionNode and
 * ~GenNode. Serial tree nodes are not using the destructors, but explicitely call to deallocNodes or deallocGenNodes
 * Gen nodes and loose nodes are not counted with MWTree->[in/de]crementNodeCount()
 */
template <int D>
GenNodeAllocator<D>::GenNodeAllocator(FunctionTree<D> *tree)
        : NodeAllocator<D>(tree, nullptr)
        , lastNode(nullptr) {
    this->coeffsPerNode = this->tree_p->getKp1_d();
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
    this->lastNode = static_cast<FunctionNode<D> *>(this->sNodes);

    MRCPP_INIT_OMP_LOCK();
}

template <int D> GenNodeAllocator<D>::~GenNodeAllocator() {
    for (int i = 0; i < this->nodeChunks.size(); i++) delete[](char *)(this->nodeChunks[i]);
    if (not this->isShared()) // if the data is shared, it must be freed by MPI_Win_free
        for (int i = 0; i < this->nodeCoeffChunks.size(); i++) delete[] this->nodeCoeffChunks[i];

    this->nodeStackStatus.clear();

    MRCPP_DESTROY_OMP_LOCK();
}

template <int D> void GenNodeAllocator<D>::allocRoots(MWTree<D> &tree) {
    NOT_REACHED_ABORT;
}

template <int D> void GenNodeAllocator<D>::allocChildren(MWNode<D> &parent, bool allocCoefs) {
    // NB: serial tree MUST generate all children consecutively
    // all children must be generated at once if several threads are active
    int sIx;
    int nChildren = parent.getTDim();
    double *coefs_p = nullptr;
    FunctionNode<D> *child_p = nullptr;
    if (allocCoefs) {
        child_p = this->allocNodes(nChildren, &sIx, &coefs_p);
    } else {
        child_p = this->allocNodes(nChildren, &sIx);
    }

    // position of first child
    parent.childSerialIx = sIx;
    for (int cIdx = 0; cIdx < nChildren; cIdx++) {
        new (child_p) FunctionNode<D>(parent, cIdx);
        parent.children[cIdx] = child_p;

        child_p->serialIx = sIx;
        child_p->parentSerialIx = parent.serialIx;
        child_p->childSerialIx = -1;
        child_p->setIsGenNode();

        if (allocCoefs) {
            child_p->n_coefs = this->coeffsPerNode;
            child_p->coefs = coefs_p;
            child_p->setIsAllocated();
        }

        sIx++;
        child_p++;
        if (allocCoefs) coefs_p += this->coeffsPerNode;
    }
}

// return pointer to the last active node or NULL if failed
template <int D> FunctionNode<D> *GenNodeAllocator<D>::allocNodes(int nAlloc, int *serialIx, double **coefs_p) {
    MRCPP_SET_OMP_LOCK();
    *serialIx = this->topStack;
    int chunkIx = *serialIx % (this->maxNodesPerChunk);

    if (chunkIx == 0 or chunkIx + nAlloc > this->maxNodesPerChunk) {
        // we want nodes allocated simultaneously to be allocated on the same piece.
        // possibly jump over the last nodes from the old chunk
        this->topStack = this->maxNodesPerChunk * ((this->topStack + nAlloc - 1) / this->maxNodesPerChunk); // start of next chunk

        int chunk = this->topStack / this->maxNodesPerChunk; // find the right chunk

        // careful: nodeChunks.size() is an unsigned int
        if (chunk + 1 > this->nodeChunks.size()) {
            // need to allocate new chunk
            double *sNodesCoeff;
            if (this->isShared()) {
                // for coefficients, take from the shared memory block
                sNodesCoeff = this->shmem_p->sh_end_ptr;
                this->shmem_p->sh_end_ptr += (this->coeffsPerNode * this->maxNodesPerChunk);
                // may increase size dynamically in the future
                if (this->shmem_p->sh_max_ptr < this->shmem_p->sh_end_ptr) MSG_ABORT("Shared block too small");
            } else {
                sNodesCoeff = new double[this->coeffsPerNode * this->maxNodesPerChunk];
            }

            this->nodeCoeffChunks.push_back(sNodesCoeff);
            this->sNodes = (FunctionNode<D> *)new char[this->maxNodesPerChunk * sizeof(FunctionNode<D>)];
            for (int i = 0; i < this->maxNodesPerChunk; i++) {
                this->sNodes[i].serialIx = -1;
                this->sNodes[i].parentSerialIx = -1;
                this->sNodes[i].childSerialIx = -1;
            }
            this->nodeChunks.push_back(this->sNodes);

            // allocate new chunk in nodeStackStatus
            int oldsize = this->nodeStackStatus.size();
            int newsize = oldsize + this->maxNodesPerChunk;
            for (int i = oldsize; i < newsize; i++) this->nodeStackStatus.push_back(0);

            if (chunk % 100 == 99 and D == 3)
                println(10,
                        std::endl
                            << " number of nodes " << this->topStack << ",number of Nodechunks now "
                            << this->nodeChunks.size() << ", total size coeff  (MB) "
                            << (this->topStack * this->coeffsPerNode) / 1024 / 128);
        }
        this->lastNode = this->nodeChunks[chunk] + this->topStack % (this->maxNodesPerChunk);
        *serialIx = this->topStack;
        chunkIx = *serialIx % (this->maxNodesPerChunk);
    }
    assert((this->topStack + nAlloc - 1) / this->maxNodesPerChunk < this->nodeChunks.size());

    FunctionNode<D> *newNode = this->lastNode;
    FunctionNode<D> *newNode_cp = newNode;

    int chunk = this->topStack / this->maxNodesPerChunk; // find the right chunk
    *coefs_p = this->nodeCoeffChunks[chunk] + chunkIx * this->coeffsPerNode;

    for (int i = 0; i < nAlloc; i++) {
        if (this->nodeStackStatus[*serialIx + i] != 0)
            println(0, *serialIx + i << " NodeStackStatus: not available " << this->nodeStackStatus[*serialIx + i]);
        this->nodeStackStatus[*serialIx + i] = 1;
        newNode_cp++;
    }
    this->nNodes += nAlloc;
    this->topStack += nAlloc;
    this->lastNode += nAlloc;

    MRCPP_UNSET_OMP_LOCK();
    return newNode;
}

// return pointer to the last active node or NULL if failed
// Will not allocate coefficients
template <int D> FunctionNode<D> *GenNodeAllocator<D>::allocNodes(int nAlloc, int *serialIx) {
    MRCPP_SET_OMP_LOCK();
    *serialIx = this->topStack;
    int chunkIx = *serialIx % (this->maxNodesPerChunk);

    if (chunkIx == 0 or chunkIx + nAlloc > this->maxNodesPerChunk) {
        // we want nodes allocated simultaneously to be allocated on the same piece.
        // possibly jump over the last nodes from the old chunk
        this->topStack = this->maxNodesPerChunk * ((this->topStack + nAlloc - 1) / this->maxNodesPerChunk); // start of next chunk

        int chunk = this->topStack / this->maxNodesPerChunk; // find the right chunk

        // careful: nodeChunks.size() is an unsigned int
        if (chunk + 1 > this->nodeChunks.size()) {
            // need to allocate new chunk
            this->sNodes = (FunctionNode<D> *)new char[this->maxNodesPerChunk * sizeof(FunctionNode<D>)];
            for (int i = 0; i < this->maxNodesPerChunk; i++) {
                this->sNodes[i].serialIx = -1;
                this->sNodes[i].parentSerialIx = -1;
                this->sNodes[i].childSerialIx = -1;
            }
            this->nodeChunks.push_back(this->sNodes);

            // allocate new chunk in nodeStackStatus
            int oldsize = this->nodeStackStatus.size();
            int newsize = oldsize + this->maxNodesPerChunk;
            for (int i = oldsize; i < newsize; i++) this->nodeStackStatus.push_back(0);

            if (chunk % 100 == 99 and D == 3)
                println(10,
                        std::endl
                            << " number of nodes " << this->topStack << ",number of Nodechunks now "
                            << this->nodeChunks.size() << ", total size coeff  (MB) "
                            << (this->topStack * this->coeffsPerNode) / 1024 / 128);
        }
        this->lastNode = this->nodeChunks[chunk] + this->topStack % (this->maxNodesPerChunk);
        *serialIx = this->topStack;
        chunkIx = *serialIx % (this->maxNodesPerChunk);
    }
    assert((this->topStack + nAlloc - 1) / this->maxNodesPerChunk < this->nodeChunks.size());

    FunctionNode<D> *newNode = this->lastNode;
    FunctionNode<D> *newNode_cp = newNode;

    int chunk = this->topStack / this->maxNodesPerChunk; // find the right chunk

    for (int i = 0; i < nAlloc; i++) {
        if (this->nodeStackStatus[*serialIx + i] != 0)
            println(0, *serialIx + i << " NodeStackStatus: not available " << this->nodeStackStatus[*serialIx + i]);
        this->nodeStackStatus[*serialIx + i] = 1;
        newNode_cp++;
    }
    this->nNodes += nAlloc;
    this->topStack += nAlloc;
    this->lastNode += nAlloc;

    MRCPP_UNSET_OMP_LOCK();
    return newNode;
}

template <int D> void GenNodeAllocator<D>::deallocNodes(int serialIx) {
    MRCPP_SET_OMP_LOCK();
    this->nodeStackStatus[serialIx] = 0; // mark as available
    if (serialIx == this->topStack - 1) {  // top of stack
        while (this->nodeStackStatus[this->topStack - 1] == 0) {
            this->topStack--;
            if (this->topStack < 1) break;
        }
        // has to redefine lastNode
        int chunk = this->topStack / this->maxNodesPerChunk; // find the right chunk
        this->lastNode = this->nodeChunks[chunk] + this->topStack % (this->maxNodesPerChunk);
    }
    this->nNodes--;
    MRCPP_UNSET_OMP_LOCK();
}

template class GenNodeAllocator<1>;
template class GenNodeAllocator<2>;
template class GenNodeAllocator<3>;

} // namespace mrcpp
