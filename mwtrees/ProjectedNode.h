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
    friend class MWNode<D>;
    friend class FunctionTree<D>;
    friend class TreeAllocator<D>;

protected:
    ProjectedNode(FunctionTree<D> &t, const NodeIndex<D> &nIdx);
    ProjectedNode(ProjectedNode<D> &p, int cIdx);
    ProjectedNode(const ProjectedNode<D> &n) : FunctionNode<D>(n) { NOT_IMPLEMENTED_ABORT; }
    ProjectedNode& operator=(const ProjectedNode<D> &n) { NOT_IMPLEMENTED_ABORT; }
    virtual ~ProjectedNode();

    void createChildren() { MWNode<D>::createChildren(); this->clearIsEndNode(); }
    void deleteChildren() { MWNode<D>::deleteChildren(); this->setIsEndNode(); }

private:
    void createChild(int i);
    void genChild(int i);

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & boost::serialization::base_object<FunctionNode<D> >(*this);
    }
};

#endif /* PROJECTEDNODE_H_ */
