
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

    friend class SerialTree<D>;

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
};

#endif /* GENNODE_H_ */
