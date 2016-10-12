/**
 *
 *  \date Aug 14, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 *
 */

#ifndef PROJECTEDNODE_H_
#define PROJECTEDNODE_H_

#include "FunctionNode.h"

template<int D>
class ProjectedNode: public FunctionNode<D> {
public:
    friend class SerialTree<D>;

protected:
    ProjectedNode() : FunctionNode<D>() { }
    virtual ~ProjectedNode() { assert(this->tree == 0); }

    void createChildren() { MWNode<D>::createChildren(); this->clearIsEndNode(); }
    void deleteChildren() { MWNode<D>::deleteChildren(); this->setIsEndNode(); }

    void reCompress();
};

#endif /* PROJECTEDNODE_H_ */
