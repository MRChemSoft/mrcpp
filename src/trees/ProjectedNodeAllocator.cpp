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

#include "ProjectedNodeAllocator.h"
#include "FunctionTree.h"
#include "ProjectedNode.h"
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
ProjectedNodeAllocator<D>::ProjectedNodeAllocator(FunctionTree<D> *tree, SharedMemory *mem)
        : NodeAllocator<D>(tree, mem)
        , lastNode(nullptr) {

    this->coeffsPerNode = (1 << D) * this->tree_p->getKp1_d(); // TDim  blocks
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
    this->lastNode = static_cast<ProjectedNode<D> *>(this->sNodes);

    // make virtual table pointers
    auto *tmpNode = new ProjectedNode<D>();
    this->cvptr_ProjectedNode = *(char **)(tmpNode);
    delete tmpNode;

    MRCPP_INIT_OMP_LOCK();
}

template <int D> ProjectedNodeAllocator<D>::~ProjectedNodeAllocator() {
    for (int i = 0; i < this->nodeChunks.size(); i++) delete[](char *)(this->nodeChunks[i]);
    if (not this->isShared()) // if the data is shared, it must be freed by MPI_Win_free
        for (int i = 0; i < this->nodeCoeffChunks.size(); i++) delete[] this->nodeCoeffChunks[i];

    this->nodeStackStatus.clear();

    MRCPP_DESTROY_OMP_LOCK();
}

/** reset the start node counter */
template <int D> void ProjectedNodeAllocator<D>::clear(int n) {
    for (int i = n; i < this->nodeStackStatus.size(); i++) this->nodeStackStatus[i] = 0;
    this->nNodes = n;
    int chunk = n / this->maxNodesPerChunk;
    this->lastNode = this->nodeChunks[chunk] + n % (this->maxNodesPerChunk);

    if (this->isShared()) {
        this->shMem->sh_end_ptr =
            this->shMem->sh_start_ptr + (n / this->maxNodesPerChunk + 1) * this->coeffsPerNode * this->maxNodesPerChunk;
    }
}

template <int D> void ProjectedNodeAllocator<D>::allocRoots(MWTree<D> &tree) {
    int sIx;
    double *coefs_p;
    // reserve place for nRoots
    int nRoots = tree.getRootBox().size();
    ProjectedNode<D> *root_p = this->allocNodes(nRoots, &sIx, &coefs_p);

    MWNode<D> **roots = tree.getRootBox().getNodes();
    for (int rIdx = 0; rIdx < nRoots; rIdx++) {
        roots[rIdx] = root_p;

        *(char **)(root_p) = this->cvptr_ProjectedNode;

        root_p->tree = &tree;
        root_p->parent = nullptr;
        for (int i = 0; i < root_p->getTDim(); i++) { root_p->children[i] = nullptr; }

        root_p->maxSquareNorm = -1.0;
        root_p->maxWSquareNorm = -1.0;

        root_p->nodeIndex = tree.getRootBox().getNodeIndex(rIdx);
        root_p->hilbertPath = HilbertPath<D>();

        root_p->n_coefs = this->coeffsPerNode;
        root_p->coefs = coefs_p;

        root_p->lockX = 0;
        root_p->serialIx = sIx;
        root_p->parentSerialIx = -1; // to indicate rootnode
        root_p->childSerialIx = -1;

        root_p->status = 0;

        root_p->clearNorms();
        root_p->setIsLeafNode();
        root_p->setIsAllocated();
        root_p->clearHasCoefs();
        root_p->setIsEndNode();
        root_p->setIsRootNode();

        tree.incrementNodeCount(root_p->getScale());

        sIx++;
        root_p++;
        coefs_p += this->coeffsPerNode;
    }
}

template <int D> void ProjectedNodeAllocator<D>::allocChildren(MWNode<D> &parent) {
    int sIx;
    double *coefs_p;
    // NB: serial tree MUST generate all children consecutively
    // all children must be generated at once if several threads are active
    int nChildren = parent.getTDim();
    ProjectedNode<D> *child_p = this->allocNodes(nChildren, &sIx, &coefs_p);

    // position of first child
    parent.childSerialIx = sIx;
    for (int cIdx = 0; cIdx < nChildren; cIdx++) {
        parent.children[cIdx] = child_p;

        *(char **)(child_p) = this->cvptr_ProjectedNode;

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
        child_p->setIsEndNode();

        child_p->tree->incrementNodeCount(child_p->getScale());

        sIx++;
        child_p++;
        coefs_p += this->coeffsPerNode;
    }
}

template <int D> void ProjectedNodeAllocator<D>::allocChildrenNoCoeff(MWNode<D> &parent) {
    int sIx;
    // all children must be generated at once if several threads are active
    int nChildren = parent.getTDim();
    ProjectedNode<D> *child_p = this->allocNodes(nChildren, &sIx);

    // position of first child
    parent.childSerialIx = sIx;
    for (int cIdx = 0; cIdx < nChildren; cIdx++) {
        parent.children[cIdx] = child_p;

        *(char **)(child_p) = this->cvptr_ProjectedNode;

        child_p->tree = parent.tree;
        child_p->parent = &parent;
        for (int i = 0; i < child_p->getTDim(); i++) { child_p->children[i] = nullptr; }

        child_p->maxSquareNorm = -1.0;
        child_p->maxWSquareNorm = -1.0;

        child_p->nodeIndex = parent.getNodeIndex().child(cIdx);
        child_p->hilbertPath = HilbertPath<D>(parent.getHilbertPath(), cIdx);

        child_p->n_coefs = 0;
        child_p->coefs = nullptr;

        child_p->lockX = 0;
        child_p->serialIx = sIx;
        child_p->parentSerialIx = parent.serialIx;
        child_p->childSerialIx = -1;

        child_p->status = 0;

        child_p->clearNorms();
        child_p->setIsLeafNode();
        child_p->setIsAllocated();
        child_p->clearHasCoefs();
        child_p->setIsEndNode();

        child_p->tree->incrementNodeCount(child_p->getScale());

        sIx++;
        child_p++;
    }
}

// return pointer to the last active node or NULL if failed
template <int D> ProjectedNode<D> *ProjectedNodeAllocator<D>::allocNodes(int nAlloc, int *serialIx, double **coefs_p) {
    *serialIx = this->nNodes;
    int chunkIx = *serialIx % (this->maxNodesPerChunk);

    if (chunkIx == 0 or chunkIx + nAlloc > this->maxNodesPerChunk) {
        // we want nodes allocated simultaneously to be allocated on the same piece.
        // possibly jump over the last nodes from the old chunk
        this->nNodes = this->maxNodesPerChunk * ((this->nNodes + nAlloc - 1) / this->maxNodesPerChunk); // start of next chunk

        int chunk = this->nNodes / this->maxNodesPerChunk; // find the right chunk

        // careful: nodeChunks.size() is an unsigned int
        if (chunk + 1 > this->nodeChunks.size()) {
            // need to allocate new chunk
            double *sNodesCoeff;
            if (this->isShared()) {
                // for coefficients, take from the shared memory block
                sNodesCoeff = this->shMem->sh_end_ptr;
                this->shMem->sh_end_ptr += (this->coeffsPerNode * this->maxNodesPerChunk);
                // may increase size dynamically in the future
                if (this->shMem->sh_max_ptr < this->shMem->sh_end_ptr) MSG_ABORT("Shared block too small");
            } else {
                sNodesCoeff = new double[this->coeffsPerNode * this->maxNodesPerChunk];
            }

            this->nodeCoeffChunks.push_back(sNodesCoeff);
            this->sNodes = (ProjectedNode<D> *)new char[this->maxNodesPerChunk * sizeof(ProjectedNode<D>)];
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
            this->maxNodes = newsize;

            if (chunk % 100 == 99 and D == 3)
                println(10,
                        std::endl
                            << " number of nodes " << this->nNodes << ",number of Nodechunks now "
                            << this->nodeChunks.size() << ", total size coeff  (MB) "
                            << (this->nNodes * this->coeffsPerNode) / 1024 / 128);
        }
        this->lastNode = this->nodeChunks[chunk] + this->nNodes % (this->maxNodesPerChunk);
        *serialIx = this->nNodes;
        chunkIx = *serialIx % (this->maxNodesPerChunk);
    }
    assert((this->nNodes + nAlloc - 1) / this->maxNodesPerChunk < this->nodeChunks.size());

    ProjectedNode<D> *newNode = this->lastNode;
    ProjectedNode<D> *newNode_cp = newNode;

    int chunk = this->nNodes / this->maxNodesPerChunk; // find the right chunk
    *coefs_p = this->nodeCoeffChunks[chunk] + chunkIx * this->coeffsPerNode;

    for (int i = 0; i < nAlloc; i++) {
        if (this->nodeStackStatus[*serialIx + i] != 0)
            println(0, *serialIx + i << " NodeStackStatus: not available " << this->nodeStackStatus[*serialIx + i]);
        this->nodeStackStatus[*serialIx + i] = 1;
        newNode_cp++;
    }
    this->nNodes += nAlloc;
    this->lastNode += nAlloc;

    return newNode;
}

// return pointer to the last active node or NULL if failed
// Will not allocate coefficients
template <int D> ProjectedNode<D> *ProjectedNodeAllocator<D>::allocNodes(int nAlloc, int *serialIx) {
    *serialIx = this->nNodes;
    int chunkIx = *serialIx % (this->maxNodesPerChunk);

    if (chunkIx == 0 or chunkIx + nAlloc > this->maxNodesPerChunk) {
        // we want nodes allocated simultaneously to be allocated on the same piece.
        // possibly jump over the last nodes from the old chunk
        this->nNodes = this->maxNodesPerChunk * ((this->nNodes + nAlloc - 1) / this->maxNodesPerChunk); // start of next chunk

        int chunk = this->nNodes / this->maxNodesPerChunk; // find the right chunk

        // careful: nodeChunks.size() is an unsigned int
        if (chunk + 1 > this->nodeChunks.size()) {
            // need to allocate new chunk
            this->sNodes = (ProjectedNode<D> *)new char[this->maxNodesPerChunk * sizeof(ProjectedNode<D>)];
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
            this->maxNodes = newsize;

            if (chunk % 100 == 99 and D == 3)
                println(10,
                        std::endl
                            << " number of nodes " << this->nNodes << ",number of Nodechunks now "
                            << this->nodeChunks.size() << ", total size coeff  (MB) "
                            << (this->nNodes * this->coeffsPerNode) / 1024 / 128);
        }
        this->lastNode = this->nodeChunks[chunk] + this->nNodes % (this->maxNodesPerChunk);
        *serialIx = this->nNodes;
        chunkIx = *serialIx % (this->maxNodesPerChunk);
    }
    assert((this->nNodes + nAlloc - 1) / this->maxNodesPerChunk < this->nodeChunks.size());

    ProjectedNode<D> *newNode = this->lastNode;
    ProjectedNode<D> *newNode_cp = newNode;

    int chunk = this->nNodes / this->maxNodesPerChunk; // find the right chunk

    for (int i = 0; i < nAlloc; i++) {
        if (this->nodeStackStatus[*serialIx + i] != 0)
            println(0, *serialIx + i << " NodeStackStatus: not available " << this->nodeStackStatus[*serialIx + i]);
        this->nodeStackStatus[*serialIx + i] = 1;
        newNode_cp++;
    }
    this->nNodes += nAlloc;
    this->lastNode += nAlloc;

    return newNode;
}

template <int D> void ProjectedNodeAllocator<D>::deallocNodes(int serialIx) {
    if (this->nNodes < 0) {
        println(0, "minNodes exceeded " << this->nNodes);
        this->nNodes++;
    }
    this->nodeStackStatus[serialIx] = 0; // mark as available
    if (serialIx == this->nNodes - 1) {  // top of stack
        int topStack = this->nNodes;
        while (this->nodeStackStatus[topStack - 1] == 0) {
            topStack--;
            if (topStack < 1) break;
        }
        this->nNodes = topStack; // move top of stack
        // has to redefine lastNode
        int chunk = this->nNodes / this->maxNodesPerChunk; // find the right chunk
        this->lastNode = this->nodeChunks[chunk] + this->nNodes % (this->maxNodesPerChunk);
    }
}

/** Fill all holes in the chunks with occupied nodes, then remove all empty chunks */
template <int D> int ProjectedNodeAllocator<D>::shrinkChunks() {
    int nAlloc = (1 << D);
    if (this->maxNodesPerChunk * this->nodeChunks.size() <=
        this->getTree()->getNNodes() + this->maxNodesPerChunk + nAlloc - 1) {
        return 0; // no chunks to remove
    }

    double *coefs_p;
    MWTree<D> *tree = this->tree_p;
    int posavail = tree->getRootBox().size(); // start after root nodes
    int posocc = 0;
    int nChunksStart = this->nodeChunks.size();

    while (true) {
        // find first available spot
        // Note that last positions on a chunk cannot be used if there is no place for 8 siblings on the same chunk
        while ((this->nodeStackStatus[posavail] != 0 or
                (posavail + nAlloc - 1) / this->maxNodesPerChunk != posavail / this->maxNodesPerChunk) and
               posavail < this->nNodes)
            posavail++;
        if (posavail >= this->nNodes) break; // treated all nodes

        // find next allocated spot after available
        posocc = posavail;
        while (this->nodeStackStatus[posocc] == 0 and posavail < this->nNodes) posocc++;
        if (posocc >= this->nNodes) break; // treated all nodes

        // move node from posocc to posavail
        int ichunk = posocc / this->maxNodesPerChunk;
        int inode = posocc % this->maxNodesPerChunk;
        ProjectedNode<D> *NodeOcc = this->nodeChunks[ichunk] + inode;

        // check that all siblings are consecutive. Should never be root node.
        for (int i = 0; i < nAlloc; i++) assert(this->nodeStackStatus[posavail + i] == 0);
        for (int i = 1; i < nAlloc; i++)
            assert((NodeOcc + i)->parent->serialIx == (NodeOcc)->parent->serialIx); // siblings

        ichunk = posavail / this->maxNodesPerChunk;
        inode = posavail % this->maxNodesPerChunk;
        ProjectedNode<D> *NodeAvail = this->nodeChunks[ichunk] + inode;

        // just copy everything "as is"
        for (int i = 0; i < nAlloc * sizeof(ProjectedNode<D>); i++) ((char *)NodeAvail)[i] = ((char *)NodeOcc)[i];

        // coefs have new adresses
        for (int i = 0; i < nAlloc; i++)
            (NodeAvail + i)->coefs = this->nodeCoeffChunks[ichunk] + (inode + i) * this->coeffsPerNode;
        // copy coefs to new adress
        if (not this->isShared()) {
            for (int i = 0; i < nAlloc * this->coeffsPerNode; i++) NodeAvail->coefs[i] = NodeOcc->coefs[i];
        } else {
            if (this->shMem->rank == 0) // only master copy the data. careful with sync
                for (int i = 0; i < nAlloc * this->coeffsPerNode; i++) NodeAvail->coefs[i] = NodeOcc->coefs[i];
        }

        // new nodes have another adress
        for (int i = 0; i < nAlloc; i++) (NodeAvail + i)->serialIx = posavail + i;

        // new nodes have another adress. Update parent
        NodeAvail->parent->childSerialIx = posavail;
        for (int i = 0; i < nAlloc; i++) NodeAvail->parent->children[i] = NodeAvail + i;

        // Update children too
        for (int i = 0; i < nAlloc; i++) {
            for (int j = 0; j < (NodeAvail + i)->getNChildren(); j++)
                (NodeAvail + i)->children[j]->parentSerialIx = posavail + i;
            for (int j = 0; j < (NodeAvail + i)->getNChildren(); j++)
                (NodeAvail + i)->children[j]->parent = NodeAvail + i;
        }

        // mark moved nodes as occupied
        for (int i = 0; i < nAlloc; i++) this->nodeStackStatus[posavail + i] = 1;
        posavail += nAlloc;

        // delete "old" nodes
        for (int i = 0; i < nAlloc; i++) this->nodeStackStatus[posocc + i] = 0;
        for (int i = 0; i < nAlloc; i++) (NodeOcc + i)->serialIx = -1;
    }

    // find the last used node
    posocc = this->nNodes - 1;
    while (this->nodeStackStatus[posocc] == 0 and posocc > 0) posocc--;
    this->nNodes = posocc + 1;
    int ichunk = this->nNodes / this->maxNodesPerChunk;
    int inode = this->nNodes % this->maxNodesPerChunk;
    this->lastNode = this->nodeChunks[ichunk] + inode;

    int nChunks = posocc / this->maxNodesPerChunk + 1; // number of occupied chunks
    for (int i = nChunks; i < this->nodeChunks.size(); i++)
        delete[](char *)(this->nodeChunks[i]); // remove unused chunks

    if (this->isShared()) {
        // shared coefficients cannot be fully deallocated, only pointer is moved.
        this->shMem->sh_end_ptr -= (nChunksStart - nChunks) * this->coeffsPerNode * this->maxNodesPerChunk;
    } else {
        for (int i = nChunks; i < this->nodeCoeffChunks.size(); i++) delete[] this->nodeCoeffChunks[i];
    }

    // shrink the stacks
    this->nodeChunks.resize(nChunks);
    this->nodeCoeffChunks.resize(nChunks);
    this->nodeStackStatus.resize(nChunks * this->maxNodesPerChunk);

    this->maxNodes = this->nodeStackStatus.size();

    this->getTree()->resetEndNodeTable();

    return nChunksStart - nChunks;
}

/** Traverse tree and redefine pointer, counter and tables. */
template <int D> void ProjectedNodeAllocator<D>::rewritePointers(bool coeff) {
    this->getTree()->nNodes = 0;
    this->getTree()->nodesAtDepth.clear();
    this->getTree()->squareNorm = 0.0;
    this->getTree()->clearEndNodeTable();

    // reinitialize stacks
    for (int i = 0; i < this->nodeStackStatus.size(); i++) this->nodeStackStatus[i] = 0;
    // nodeChunks have been adapted to receiving tree and maybe larger than nodeStackStatus
    int nodecount = this->nodeChunks.size() * this->maxNodesPerChunk;
    while (nodecount > this->nodeStackStatus.size()) this->nodeStackStatus.push_back(0);
    this->maxNodes = this->nodeStackStatus.size();

    NodeBox<D> &rBox = this->getTree()->getRootBox();
    MWNode<D> **roots = rBox.getNodes();

    int DepthMax = 100, slen = 0;
    ProjectedNode<D> *stack[DepthMax * 8];
    for (int rIdx = 0; rIdx < rBox.size(); rIdx++) {
        roots[rIdx] = (this->nodeChunks[0]) + rIdx; // adress of roots are at start of NodeChunks[0] array
        stack[slen++] = (this->nodeChunks[0]) + rIdx;
    }
    this->nNodes = 0;
    while (slen) {
        ProjectedNode<D> *node = stack[--slen];
        for (int i = 0; i < node->getNChildren(); i++) {
            int n_ichunk = (node->childSerialIx + i) / this->maxNodesPerChunk;
            int n_inode = (node->childSerialIx + i) % this->maxNodesPerChunk;
            stack[slen++] = this->nodeChunks[n_ichunk] + n_inode;
        }
        int ichunk = node->serialIx / this->maxNodesPerChunk;
        int inode = node->serialIx % this->maxNodesPerChunk;

        this->nNodes = std::max(this->nNodes, ichunk * this->maxNodesPerChunk + inode + 1);
        this->getTree()->incrementNodeCount(node->getScale());
        if (node->isEndNode()) this->getTree()->squareNorm += node->getSquareNorm();
        if (node->isEndNode()) this->getTree()->endNodeTable.push_back(node);

        // normally (intel) the virtual table does not change, but we overwrite anyway
        *(char **)(node) = this->cvptr_ProjectedNode;

        node->tree = this->getTree();

        //"adress" of coefs is the same as node, but in another array
        node->coefs = nullptr;
        if (coeff) node->coefs = this->nodeCoeffChunks[ichunk] + inode * this->coeffsPerNode;

        // adress of parent and children must be corrected
        // can be on a different chunks
        if (node->parentSerialIx >= 0) {
            int n_ichunk = node->parentSerialIx / this->maxNodesPerChunk;
            int n_inode = node->parentSerialIx % this->maxNodesPerChunk;
            node->parent = this->nodeChunks[n_ichunk] + n_inode;
            assert(node->parent->serialIx == node->parentSerialIx);
        } else {
            node->parent = nullptr;
        }

        for (int i = 0; i < node->getNChildren(); i++) {
            int n_ichunk = (node->childSerialIx + i) / this->maxNodesPerChunk;
            int n_inode = (node->childSerialIx + i) % this->maxNodesPerChunk;
            node->children[i] = this->nodeChunks[n_ichunk] + n_inode;
        }
        this->nodeStackStatus[node->serialIx] = 1; // occupied
    }
    int ichunk = this->nNodes / this->maxNodesPerChunk;
    int inode = this->nNodes % this->maxNodesPerChunk;
    this->lastNode = this->nodeChunks[ichunk] + inode;
}

template <int D> int ProjectedNodeAllocator<D>::getNChunksUsed() const {
    return (this->nNodes + this->maxNodesPerChunk - 1) / this->maxNodesPerChunk;
}

template <int D> void ProjectedNodeAllocator<D>::print() const {
    int n = 0;
    for (int iChunk = 0; iChunk < getNChunks(); iChunk++) {
        int iShift = iChunk * this->maxNodesPerChunk;
        printout(0, "new chunk \n");
        for (int i = 0; i < this->maxNodesPerChunk; i++) {
            int status = this->nodeStackStatus[iShift + i];
            int sIdx = this->nodeChunks[iChunk][i].serialIx;
            int pIdx = this->nodeChunks[iChunk][i].parentSerialIx;
            int cIdx = this->nodeChunks[iChunk][i].childSerialIx;
            printout(0, std::setw(4) << n++);
            printout(0, std::setw(4) << status);
            printout(0, std::setw(6) << sIdx);
            printout(0, std::setw(6) << pIdx);
            printout(0, std::setw(6) << cIdx << "   ");
            if (status == 1) printout(0, this->nodeChunks[iChunk][i].getNodeIndex());
            printout(0, "\n");
        }
    }
}

template <int D> void ProjectedNodeAllocator<D>::initChunk(int iChunk, bool coeff) {
    if (iChunk < getNChunks()) {
        this->sNodes = this->nodeChunks[iChunk];
    } else {
        double *sNodesCoeff = nullptr;
        if (this->isShared()) {
            if (iChunk == 0) this->shMem->sh_end_ptr = this->shMem->sh_start_ptr;
            sNodesCoeff = this->shMem->sh_end_ptr;
            this->shMem->sh_end_ptr += (this->maxNodesPerChunk * this->coeffsPerNode);
            if (this->shMem->sh_end_ptr > this->shMem->sh_max_ptr) MSG_ABORT("Shared block too small");
        } else if (coeff) {
            sNodesCoeff = new double[getCoeffChunkSize()];
        }
        if (coeff) this->nodeCoeffChunks.push_back(sNodesCoeff);
        this->sNodes = (ProjectedNode<D> *)new char[getNodeChunkSize()];
        this->nodeChunks.push_back(this->sNodes);
    }
}

template class ProjectedNodeAllocator<1>;
template class ProjectedNodeAllocator<2>;
template class ProjectedNodeAllocator<3>;

} // namespace mrcpp
