/**
 *
 *  \date July, 2016
 *  \author Peter Wind <peter.wind@uit.no> \n
 *  CTCC, University of Tromsø
 *
 */

#pragma once

#pragma GCC system_header
#include <Eigen/Core>
#include <vector>

#include "mrcpp_declarations.h"

namespace mrcpp {

template<int D>
class SerialTree {
public:
    SerialTree(MWTree<D> *tree, SharedMemory *mem);
    virtual ~SerialTree() { }

    MWTree<D>* getTree() { return this->tree_p; }
    SharedMemory* getMemory() { return this->shMem; }

    bool isShared() const { if (this->shMem == 0) return false; return true; }

    virtual void allocRoots(MWTree<D> &tree) = 0;
    virtual void allocChildren(MWNode<D> &parent) = 0;
    virtual void allocGenChildren(MWNode<D> &parent) = 0;

    virtual void deallocNodes(int serialIx) = 0;
    virtual void deallocGenNodes(int serialIx) = 0;
    virtual void deallocGenNodeChunks() = 0;

    void S_mwTransform(double* coeff_in, double* coeff_out, bool readOnlyScaling, int stride, bool overwrite = true);
    void S_mwTransformBack(double* coeff_in, double* coeff_out, int stride);

    int nNodes;                 //number of Nodes already defined
    int maxNodesPerChunk;
    std::vector<int> nodeStackStatus;
    int sizeNodeCoeff;          //size of coeff for one node
    double **coeffStack;
    int maxNodes;               //max number of nodes that can be defined

protected:
    MWTree<D> *tree_p;
    SharedMemory *shMem;
};

}
