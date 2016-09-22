/**
 *  Simple n-dimensional node
 *
 *  Created on: May 29, 2009
 *      Author: jonas
 */

#ifndef MWNODE_H_
#define MWNODE_H_

#include <Eigen/Core>

#include "parallel.h"
#include "macros.h"
#include "mrcpp_declarations.h"

#include "MWTree.h"
#include "HilbertPath.h"
#include "SerialTree.h"

#ifdef OPENMP
#define SET_NODE_LOCK() omp_set_lock(&this->node_lock)
#define UNSET_NODE_LOCK() omp_unset_lock(&this->node_lock)
#define TEST_NODE_LOCK() omp_test_lock(&this->node_lock)
#else
#define SET_NODE_LOCK()
#define UNSET_NODE_LOCK()
#define TEST_NODE_LOCK() false
#endif

template<int D>
class MWNode {
public:
    int getKp1() const { return getMWTree().getKp1(); }
    int getKp1_d() const { return getMWTree().getKp1_d(); }
    int getOrder() const { return getMWTree().getOrder(); }
    int getScalingType() const { return getMWTree().getMRA().getScalingBasis().getScalingType(); }
    int getTDim() const { return getMWTree().getTDim(); }
    int getDepth() const { return getNodeIndex().getScale()-getMWTree().getRootScale(); }
    int getScale() const { return getNodeIndex().getScale(); }
    int getNChildren() const { if (isBranchNode()) return getTDim(); return 0; }
    int getRank() const { return this->NodeRank; }
    const int *getTranslation() const { return getNodeIndex().getTranslation(); }

    const NodeIndex<D> &getNodeIndex() const { return this->nodeIndex; }
    const HilbertPath<D> &getHilbertPath() const { return this->hilbertPath; }

    void getCenter(double *r) const;
    void getBounds(double *lb, double *ub) const;

    bool hasCoord(const double *r) const;
    bool isCompatible(const MWNode<D> &node);
    bool isAncestor(const NodeIndex<D> &idx) const;
    bool isDecendant(const NodeIndex<D> &idx) const;

    inline bool hasCoefs() const;
    inline bool isRootNode() const;
    inline bool isEndNode() const;
    inline bool isGenNode() const;
    inline bool isLeafNode() const;
    inline bool isAllocated() const;
    inline bool isBranchNode() const;
    inline bool isLooseNode() const;
    inline bool hasWCoefs() const;

    double getSquareNorm() const { return this->squareNorm; }
    double getScalingNorm() const;
    virtual double getWaveletNorm() const;
    double getComponentNorm(int i) const { return this->componentNorms[i]; }
    bool hasComponentNorms() const;

    int getNCoefs() const { return this->n_coefs; }
    void getCoefs(Eigen::VectorXd &c) const;
    void printCoefs() const;

    double* getCoefs() { return this->coefs; }
    const double* getCoefs() const { return this->coefs; }

    MWTree<D>& getMWTree() { return static_cast<MWTree<D> &>(*this->tree); }
    MWNode<D>& getMWParent() { return static_cast<MWNode<D> &>(*this->parent); }
    MWNode<D>& getMWChild(int i) { return static_cast<MWNode<D> &>(*this->children[i]); }

    const MWTree<D>& getMWTree() const { return static_cast<const MWTree<D> &>(*this->tree); }
    const MWNode<D>& getMWParent() const { return static_cast<const MWNode<D> &>(*this->parent); }
    const MWNode<D>& getMWChild(int i) const { return static_cast<const MWNode<D> &>(*this->children[i]); }

    template<int T>
    friend std::ostream& operator<<(std::ostream &o, const MWNode<T> &nd);

    friend class TreeAdaptor<D>;
    friend class TreeBuilder<D>;
    friend class DefaultCalculator<D>;
    friend class ProjectionCalculator<D>;
    friend class AdditionCalculator<D>;
    friend class MultiplicationCalculator<D>;
    friend class OperApplicationCalculator<D>;
    friend class CrossCorrelationCalculator;
    friend class DerivativeCalculator;
    friend class MWTree<D>;
    friend class SerialTree<D>;
    friend class FunctionTree<D>;
    friend class GenNode<D>;
    friend class NodeBox<D>;
    friend class ProjectedNode<D>;

protected:
    MWTree<D> *tree;
    MWNode<D> *parent;	    ///< Parent node
    MWNode<D> *children[1<<D];    ///< 2^D children

    const NodeIndex<D> nodeIndex;
    const HilbertPath<D> hilbertPath;

    int NodeRank;
    int NodeCoeffIx;

    double squareNorm;
    double componentNorms[1<<D]; ///< 2^D components

    Eigen::VectorXd oldcoefvec;
    Eigen::VectorXd *oldcoefs;

    int n_coefs;
    double *coefs;

    MWNode(MWTree<D> &t, const NodeIndex<D> &nIdx);
    MWNode(MWNode<D> &p, int cIdx);
    MWNode(const MWNode<D> &n);
    MWNode& operator=(const MWNode<D> &n) { NOT_IMPLEMENTED_ABORT; }
    virtual ~MWNode();

    void setHasCoefs() { SET_BITS(status, FlagHasCoefs | FlagAllocated); }
    void setHasWCoefs() { SET_BITS(status, FlagHasWCoefs); }
    void setIsEndNode() { SET_BITS(status, FlagEndNode); }
    void setIsGenNode() { SET_BITS(status, FlagGenNode); }
    void setIsRootNode() { SET_BITS(status, FlagRootNode); }
    void setIsLeafNode() { CLEAR_BITS(status, FlagBranchNode); }
    void setIsAllocated() { SET_BITS(status, FlagAllocated); }
    void setIsBranchNode() { SET_BITS(status, FlagBranchNode); }
    void setIsLooseNode() { SET_BITS(status, FlagLooseNode); }
    void clearHasCoefs() { CLEAR_BITS(status, FlagHasCoefs);}
    void clearHasWCoefs() { CLEAR_BITS(status, FlagHasWCoefs);}
    void clearIsEndNode() { CLEAR_BITS(status, FlagEndNode); }
    void clearIsRootNode() { CLEAR_BITS(status, FlagRootNode); }
    void clearIsAllocated() { CLEAR_BITS(status, FlagAllocated); }

    virtual void allocCoefs(int n_blocks, int block_size);
    virtual void freeCoefs();

    //virtual void setCoefs(const Eigen::VectorXd &c);
    virtual void zeroCoefs();

    void setCoefBlock(int block, int block_size, const double *c);
    void addCoefBlock(int block, int block_size, const double *c);
    void zeroCoefBlock(int block, int block_size);

    void calcNorms();
    void zeroNorms();
    void clearNorms();

    double estimateError(bool absPrec);

    virtual double calcComponentNorm(int i) const;

    virtual void cvTransform(int kind);
    virtual void orig_mwTransform(int kind);
    virtual void mwTransform(int kind);

    bool crop(double prec, NodeIndexSet *cropIdx = 0);
    void reCompress(bool overwrite);

    virtual void copyChildren(const MWNode<D> &node) { NOT_IMPLEMENTED_ABORT; }
    virtual void createChildren();
    virtual void deleteChildren();
    virtual void genChildren();

    virtual void genChild(int cIdx) { NOT_IMPLEMENTED_ABORT; }
    virtual void createChild(int cIdx) { NOT_IMPLEMENTED_ABORT; }

    virtual void giveChildrenCoefs(bool overwrite = true);
    virtual void copyCoefsFromChildren();

    int getChildIndex(const NodeIndex<D> &nIdx) const;
    int getChildIndex(const double *r) const;

    void lockNode() { SET_NODE_LOCK(); }
    void unlockNode() { UNSET_NODE_LOCK(); }
    bool testLock() { return TEST_NODE_LOCK(); }

    bool diffBranch(const MWNode<D> &rhs) const;
    inline bool checkStatus(unsigned char mask) const;

    MWNode<D> *retrieveNode(const double *r, int depth);
    MWNode<D> *retrieveNode(const NodeIndex<D> &idx);

    const MWNode<D> *retrieveNodeNoGen(const NodeIndex<D> &idx) const;
    MWNode<D> *retrieveNodeNoGen(const NodeIndex<D> &idx);

    const MWNode<D> *retrieveNodeOrEndNode(const double *r, int depth) const;
    MWNode<D> *retrieveNodeOrEndNode(const double *r, int depth);

    const MWNode<D> *retrieveNodeOrEndNode(const NodeIndex<D> &idx) const;
    MWNode<D> *retrieveNodeOrEndNode(const NodeIndex<D> &idx);

    virtual void clearGenerated();
    void deleteGenerated();

    static const unsigned char FlagBranchNode = B8(00000001);
    static const unsigned char FlagGenNode    = B8(00000010);
    static const unsigned char FlagHasCoefs   = B8(00000100);
    static const unsigned char FlagAllocated  = B8(00001000);
    static const unsigned char FlagEndNode    = B8(00010000);
    static const unsigned char FlagRootNode   = B8(00100000);
    static const unsigned char FlagLooseNode  = B8(01000000);
    static const unsigned char FlagHasWCoefs  = B8(10000000);//if Wavelet coefficients are not zero
#ifdef OPENMP
    omp_lock_t node_lock;
#endif

private:
    unsigned char status;
};

/** Allocation status of s/d-coefs is stored in the status bits for
 * serialization purposes. It's not enough to test if coefs == 0.
 */
template<int D>
bool MWNode<D>::isAllocated() const {
    if (this->status & FlagAllocated) {
        return true;
    }
    return false;
}

template<int D>
bool MWNode<D>::hasCoefs() const {
    if (this->status & FlagHasCoefs) {
        return true;
    }
    return false;
}

template<int D>
bool MWNode<D>::hasWCoefs() const {
    if (this->status & FlagHasWCoefs) {
        return true;
    }
    return false;
}

template<int D>
bool MWNode<D>::isGenNode() const {
    if (this->status & FlagGenNode) {
        return true;
    }
    return false;
}

template<int D>
bool MWNode<D>::isLeafNode() const {
    if (this->status & FlagBranchNode) {
        return false;
    }
    return true;
}

template<int D>
bool MWNode<D>::isBranchNode() const {
    if (this->status & FlagBranchNode) {
        return true;
    }
    return false;
}

template<int D>
bool MWNode<D>::isLooseNode() const {
    if (this->status & FlagLooseNode) {
        return true;
    }
    return false;
}

template<int D>
bool MWNode<D>::isEndNode() const {
    if (this->status & FlagEndNode) {
        return true;
    }
    return false;
}

template<int D>
bool MWNode<D>::isRootNode() const {
    if (this->status & FlagRootNode) {
        return true;
    }
    return false;
}

template<int D>
bool MWNode<D>::checkStatus(unsigned char mask) const {
    if (mask == (this->status & mask)) {
        return true;
    }
    return false;
}

template<int D>
std::ostream& operator<<(std::ostream &o, const MWNode<D> &nd) {
    std::string flags ="       ";
    o << nd.getNodeIndex();
    if (nd.isRootNode()) {
        flags[0] = 'R';
    }
    if (nd.isEndNode()) {
        flags[1] = 'E';
    }
    if (nd.isBranchNode()) {
        flags[2] = 'B';
    } else {
        flags[2] = 'L';
    }
    if (nd.isGenNode()) {
        flags[3] = 'G';
    } else {
        flags[3] = 'P';
    }
    if (nd.isAllocated()) {
        flags[4] = 'A';
    }
    if (nd.hasCoefs()) {
        flags[5] = 'C';
	if (nd.hasWCoefs()) {
	  flags[6] = 'W';
	}
    }
    o << " " << flags;
    o << " sqNorm=" << nd.squareNorm;
    if (nd.hasCoefs()) {
        o << " Coefs={";
        o << nd.getCoefs_d()[0] << ", " <<
             nd.getCoefs_d()[nd.getNCoefs() - 1] << "}";
    }
    return o;
}

#endif /* MWNODE_H_ */
