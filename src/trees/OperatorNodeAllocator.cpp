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

#include "OperatorNodeAllocator.h"
#include "OperatorNode.h"
#include "OperatorTree.h"
#include "utils/Printer.h"

namespace mrcpp {

int NOtrees = 0;

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
OperatorNodeAllocator::OperatorNodeAllocator(OperatorTree *tree)
        : NodeAllocator<2>(tree, nullptr)
        , sNodes(nullptr)
        , lastNode(nullptr) {

    this->topStack = 0;
    NOtrees++;

    this->coeffsPerNode = 4 * this->tree_p->getKp1_d();

    this->maxNodesPerChunk = 1024;
    this->lastNode = (OperatorNode *)this->sNodes; // position of last allocated node

    MRCPP_INIT_OMP_LOCK();
}

/** SerialTree destructor. */
OperatorNodeAllocator::~OperatorNodeAllocator() {
    for (auto &nodeChunk : this->nodeChunks) delete[](char *) nodeChunk;
    for (auto &nodeCoeffChunk : this->nodeCoeffChunks) delete[] nodeCoeffChunk;

    //    delete[] this->nodeStackStatus;
    this->nodeStackStatus.clear();

    NOtrees--;

    MRCPP_DESTROY_OMP_LOCK();
}

void OperatorNodeAllocator::allocRoots(MWTree<2> &tree) {
    int sIx;
    double *coefs_p;
    // reserve place for nRoots
    int nRoots = tree.getRootBox().size();
    OperatorNode *root_p = this->allocNodes(nRoots, &sIx, &coefs_p);

    MWNode<2> **roots = tree.getRootBox().getNodes();
    for (int rIdx = 0; rIdx < nRoots; rIdx++) {
        new (root_p) OperatorNode(tree, rIdx);
        roots[rIdx] = root_p;

        root_p->serialIx = sIx;
        root_p->parentSerialIx = -1; // to indicate rootnode
        root_p->childSerialIx = -1;

        root_p->n_coefs = this->coeffsPerNode;
        root_p->coefs = coefs_p;
        root_p->setIsAllocated();
        root_p->tree->incrementNodeCount(root_p->getScale());

        sIx++;
        root_p++;
        coefs_p += this->coeffsPerNode;
    }
}

void OperatorNodeAllocator::allocChildren(MWNode<2> &parent, bool allocCoefs) {
    // NB: serial tree MUST generate all children consecutively
    // all children must be generated at once if several threads are active
    int sIx;
    int nChildren = parent.getTDim();
    double *coefs_p = nullptr;
    OperatorNode *child_p = nullptr;
    if (allocCoefs) {
        child_p = this->allocNodes(nChildren, &sIx, &coefs_p);
    } else {
        NOT_REACHED_ABORT;
    }

    // position of first child
    parent.childSerialIx = sIx;
    for (int cIdx = 0; cIdx < nChildren; cIdx++) {
        new (child_p) OperatorNode(parent, cIdx);
        parent.children[cIdx] = child_p;

        child_p->serialIx = sIx;
        child_p->parentSerialIx = parent.serialIx;
        child_p->childSerialIx = -1;

        child_p->n_coefs = this->coeffsPerNode;
        child_p->coefs = coefs_p;
        child_p->setIsAllocated();
        child_p->tree->incrementNodeCount(child_p->getScale());

        sIx++;
        child_p++;
        coefs_p += this->coeffsPerNode;
    }
}

// return pointer to the last active node or NULL if failed
OperatorNode *OperatorNodeAllocator::allocNodes(int nAlloc, int *serialIx, double **coefs_p) {
    *serialIx = this->topStack;
    int chunkIx = (*serialIx) % (this->maxNodesPerChunk);

    if (chunkIx == 0 or chunkIx + nAlloc > this->maxNodesPerChunk) {
        // start on new chunk
        // we want nodes allocated simultaneously to be allocated on the same piece.
        // possibly jump over the last nodes from the old chunk
        this->topStack = this->maxNodesPerChunk * ((this->topStack + nAlloc - 1) / this->maxNodesPerChunk); // start of next chunk

        int chunk = this->topStack / this->maxNodesPerChunk; // find the right chunk

        // careful: nodeChunks.size() is an unsigned int
        if (chunk + 1 > this->nodeChunks.size()) {
            // need to allocate new chunk
            this->sNodes = (OperatorNode *)new char[this->maxNodesPerChunk * sizeof(OperatorNode)];
            this->nodeChunks.push_back(this->sNodes);
            auto *sNodesCoeff = new double[this->coeffsPerNode * this->maxNodesPerChunk];
            this->nodeCoeffChunks.push_back(sNodesCoeff);
            // allocate new chunk in nodeStackStatus
            int oldsize = this->nodeStackStatus.size();
            int newsize = oldsize + this->maxNodesPerChunk;
            for (int i = oldsize; i < newsize; i++) this->nodeStackStatus.push_back(0);
        }
        this->lastNode = this->nodeChunks[chunk] + this->topStack % (this->maxNodesPerChunk);
        *serialIx = this->topStack;
        chunkIx = *serialIx % (this->maxNodesPerChunk);
    }
    assert((this->topStack + nAlloc - 1) / this->maxNodesPerChunk < this->nodeChunks.size());

    OperatorNode *newNode = this->lastNode;
    OperatorNode *newNode_cp = newNode;

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

    return newNode;
}

void OperatorNodeAllocator::deallocNodes(int serialIx) {
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
}

} // namespace mrcpp
