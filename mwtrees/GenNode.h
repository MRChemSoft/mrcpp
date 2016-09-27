
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
    void getCoefsNoLock(Eigen::VectorXd &vec);
    void getCoefs(Eigen::VectorXd &vec);

    double getWaveletNorm() const { return 0.0; }

    friend class MWNode<D>;
    friend class ProjectedNode<D>;
    friend class SerialTree<D>;

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

    virtual void allocCoefs(int n_blocks, int block_size);
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
};


#endif /* GENNODE_H_ */
