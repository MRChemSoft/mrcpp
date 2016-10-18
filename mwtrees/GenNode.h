
/*
 *
 *  \date Oct 18, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */

#ifndef GENNODE_H_
#define GENNODE_H_

#include "FunctionNode.h"

template<int D> class SerialFunctionTree;

template<int D>
class GenNode: public FunctionNode<D> {
public:
    double getWaveletNorm() const { return 0.0; }

    virtual void createChildren() { NOT_REACHED_ABORT; }
    virtual void genChildren() {
        if (this->isBranchNode()) MSG_FATAL("Node already has children");
        this->tree->getSerialTree()->allocGenChildren(*this);
        this->setIsBranchNode();
    }
    virtual void cvTransform(int kind) { NOT_IMPLEMENTED_ABORT; }
    virtual void mwTransform(int kind) { NOT_IMPLEMENTED_ABORT; }

    friend class SerialFunctionTree<D>;

protected:
    GenNode() : FunctionNode<D>() { }
    virtual ~GenNode() { assert(this->tree == 0); }

    double calcComponentNorm(int i) const {
        if (i == 0) {
            return MWNode<D>::calcComponentNorm(0);
        } else {
            return 0.0;
        }
    }

    virtual void dealloc() {
        this->tree->decrementGenNodeCount();
        this->tree->getSerialTree()->deallocGenNodes(this->getSerialIx());
    }

    void reCompress() { NOT_IMPLEMENTED_ABORT; }
};

#endif /* GENNODE_H_ */
