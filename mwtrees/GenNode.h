
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

template<int D>
class GenNode: public FunctionNode<D> {
public:
    double getWaveletNorm() const { return 0.0; }

    virtual void cvTransform(int kind) { NOT_IMPLEMENTED_ABORT; }
    virtual void mwTransform(int kind) { NOT_IMPLEMENTED_ABORT; }

    friend class SerialTree<D>;

protected:
    GenNode() : FunctionNode<D>() { }
    virtual ~GenNode() { assert(this->tree == 0); }

    virtual void dealloc() {
	this->tree->decrementGenNodeCount();
	this->tree->serialTree_p->deallocGenNodes(this->getSerialIx());
    }

    void reCompress() { NOT_IMPLEMENTED_ABORT; }

    double calcComponentNorm(int i) const {
        if (i == 0) {
            return MWNode<D>::calcComponentNorm(0);
        } else {
            return 0.0;
        }
    }
};

#endif /* GENNODE_H_ */
