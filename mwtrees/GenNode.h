
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

template<int D> class ProjectedNode;

template<int D>
class GenNode: public FunctionNode<D> {
public:
    Eigen::VectorXd& getCoefsNoLock();
    Eigen::VectorXd& getCoefs();
    const Eigen::VectorXd& getCoefs() const;

    double getWaveletNorm() const { return 0.0; }

    friend class MWNode<D>;
    friend class ProjectedNode<D>;

protected:
    GenNode(ProjectedNode<D> &p, int cIdx);
    GenNode(GenNode<D> &p, int cIdx);
    GenNode(const GenNode<D> &n) : FunctionNode<D>(n) { NOT_IMPLEMENTED_ABORT; }
    GenNode& operator=(const GenNode<D> &n) { NOT_IMPLEMENTED_ABORT; }
    virtual ~GenNode();

    double calcComponentNorm(int i) const {
        if (i == 0) {
            return MWNode<D>::calcComponentNorm(0);
        } else {
            return 0.0;
        }
    }

    virtual void allocCoefs(int nBlocks);
    virtual void freeCoefs();
    virtual void clearGenerated();

    const ProjectedNode<D> *getGenRootNode() const { return this->genRootNode; }
    ProjectedNode<D> *getGenRootNode() { return this->genRootNode; }

private:
    ProjectedNode<D> *genRootNode;

    void createChild(int i);
    void genChild(int i);

    void lockSiblings();
    void unlockSiblings();

    void regenerateCoefs();

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & boost::serialization::base_object<FunctionNode<D> >(*this);
    }
};


#endif /* GENNODE_H_ */
