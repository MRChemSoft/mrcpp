/**
 *  Simple n-dimensional node
 *
 *  Created on: May 29, 2009
 *      Author: jonas
 */

#ifndef MWNODE_H_
#define MWNODE_H_

#include <Eigen/Core>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/utility.hpp>

#include "parallel.h"
#include "macros.h"
#include "mwrepr_declarations.h"

#include "MWTree.h"
#include "HilbertPath.h"

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
    virtual ~MWNode();

    double estimateError(bool absPrec);

    int getKp1() const { return getMWTree().getKp1(); }
    int getKp1_d() const { return getMWTree().getKp1_d(); }
    int getOrder() const { return getMWTree().getOrder(); }
    int getScalingType() const { return getMWTree().getMRA().getScalingBasis().getScalingType(); }
    int getTDim() const { return getMWTree().getTDim(); }
    int getDepth() const { return getNodeIndex().getScale()-getMWTree().getRootScale(); }
    int getScale() const { return getNodeIndex().getScale(); }
    int getRankId() const { return getNodeIndex().getRankId(); }
    int getNChildren() const { if (isBranchNode()) return getTDim(); return 0; }
    const int *getTranslation() const { return getNodeIndex().getTranslation(); }

    NodeIndex<D> &getNodeIndex() { return this->nodeIndex; }
    const NodeIndex<D> &getNodeIndex() const { return this->nodeIndex; }
    const HilbertPath<D> &getHilbertPath() const { return this->hilbertPath; }

    void getCenter(double *r) const;
    void getBounds(double *lb, double *ub) const;

    inline bool hasCoefs() const;
    inline bool isRootNode() const;
    inline bool isEndNode() const;
    inline bool isGenNode() const;
    inline bool isLeafNode() const;
    inline bool isAllocated() const;
    inline bool isBranchNode() const;

    void setHasCoefs() { SET_BITS(status, FlagHasCoefs | FlagAllocated); }
    void setIsEndNode() { SET_BITS(status, FlagEndNode); }
    void setIsGenNode() { SET_BITS(status, FlagGenNode); }
    void setIsRootNode() { SET_BITS(status, FlagRootNode); }
    void setIsLeafNode() { CLEAR_BITS(status, FlagBranchNode); }
    void setIsAllocated() { SET_BITS(status, FlagAllocated); }
    void setIsBranchNode() { SET_BITS(status, FlagBranchNode); }
    void clearHasCoefs() { CLEAR_BITS(status, FlagHasCoefs);}
    void clearIsEndNode() { CLEAR_BITS(status, FlagEndNode); }
    void clearIsRootNode() { CLEAR_BITS(status, FlagRootNode); }
    void clearIsAllocated() { CLEAR_BITS(status, FlagAllocated); }

    bool hasCoord(const double *r) const;
    bool isCompatible(const MWNode<D> &node);
    bool isAncestor(const NodeIndex<D> &idx) const;
    bool isDecendant(const NodeIndex<D> &idx) const;

    int getChildIndex(const NodeIndex<D> &nIdx) const;
    int getChildIndex(const double *r) const;

    void lockNode() { SET_NODE_LOCK(); }
    void unlockNode() { UNSET_NODE_LOCK(); }
    bool testLock() { return TEST_NODE_LOCK(); }

    void setRankId(int n) { getNodeIndex().setRankId(n); }
    bool isLocal() const {
        if (this->getRankId() == getMWTree().getRankId()) {
            return true;
        }
        return false;
    }
    bool isCommon() const {
        if (this->getRankId() < 0) {
            return true;
        }
        return false;
    }
    bool isForeign() const {
        if (isLocal() or isCommon()) {
            return false;
        }
        return true;
    }

    double getSquareNorm() const { return this->squareNorm; }
    double getScalingNorm() const { return this->componentNorms[0]; }
    double getWaveletNorm() const { return calcWaveletNorm(); }
    double getComponentNorm(int i) const { return this->componentNorms[i]; }
    bool hasComponentNorms() const;

    int getNCoefs() const { return this->coefs->size(); }
    virtual Eigen::VectorXd &getCoefs() { if (not this->isAllocated()) allocCoefs(); return *this->coefs; }
    virtual const Eigen::VectorXd &getCoefs() const { return *this->coefs; }

    virtual void setCoefs(const Eigen::VectorXd &c);
    virtual void zeroCoefs();

    virtual void cvTransform(int kind);
    virtual void mwTransform(int kind);

    MWTree<D>& getMWTree() { return static_cast<MWTree<D> &>(*this->tree); }
    MWNode<D>& getMWParent() { return static_cast<MWNode<D> &>(*this->parent); }
    MWNode<D>& getMWChild(int i) { return static_cast<MWNode<D> &>(*this->children[i]); }

    const MWTree<D>& getMWTree() const { return static_cast<const MWTree<D> &>(*this->tree); }
    const MWNode<D>& getMWParent() const { return static_cast<const MWNode<D> &>(*this->parent); }
    const MWNode<D>& getMWChild(int i) const { return static_cast<const MWNode<D> &>(*this->children[i]); }

    template<int T>
    friend std::ostream& operator<<(std::ostream &o, const MWNode<T> &nd);

    friend class TreeCalculator<D>;
    friend class ProjectionCalculator<D>;
    friend class AdditionCalculator<D>;
    friend class MWTree<D>;

protected:
    MWTree<D> *tree;
    MWNode<D> *parent;	    ///< Parent node
    MWNode<D> *children[1<<D];    ///< 2^D children

    NodeIndex<D> nodeIndex;
    const HilbertPath<D> hilbertPath;

    double squareNorm;
    double componentNorms[1<<D]; ///< 2^D components
    Eigen::VectorXd *coefs;

    MWNode(MWTree<D> &t, const NodeIndex<D> &nIdx);
    MWNode(MWNode<D> &p, int cIdx);
    MWNode(const MWNode<D> &n);
    MWNode& operator=(const MWNode<D> &n) { NOT_IMPLEMENTED_ABORT; }

    virtual void allocCoefs(int nCoefs = -1);
    virtual void freeCoefs();

    void calcNorms();
    void zeroNorms();
    void clearNorms();

    double calcSquareNorm() const;
    virtual double calcWaveletNorm() const;
    virtual double calcComponentNorm(int i) const;

    bool crop(double prec, NodeIndexSet *cropIdx = 0);
    void reCompress(bool overwrite = true);

    virtual void copyChildren(const MWNode<D> &node) { NOT_IMPLEMENTED_ABORT; }
    virtual void createChildren();
    virtual void deleteChildren();
    void genChildren();

    virtual void genChild(int cIdx) = 0;
    virtual void createChild(int cIdx) = 0;

    virtual void giveChildrenCoefs(bool overwrite = true);
    virtual void copyCoefsFromChildren(Eigen::VectorXd &c);

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

    void assignDecendantTags(int rank);
    void broadcastCoefs(int src, mpi::communicator *comm  = 0);
    virtual mpi::request isendCoefs(int who, int tag, int comp = -1);
    virtual mpi::request ireceiveCoefs(int who, int tag, int comp = -1);

    static const unsigned char FlagBranchNode = B8(00000001);
    static const unsigned char FlagGenNode    = B8(00000010);
    static const unsigned char FlagHasCoefs   = B8(00000100);
    static const unsigned char FlagAllocated  = B8(00001000);
    static const unsigned char FlagEndNode    = B8(00010000);
    static const unsigned char FlagRootNode   = B8(00100000);
#ifdef OPENMP
    omp_lock_t node_lock;
#endif

private:
    unsigned char status;

    friend class boost::serialization::access;
    template<class Archive>
    void save(Archive & ar, const unsigned int version) const {
        NOT_IMPLEMENTED_ABORT
    }
    template<class Archive>
    void load(Archive & ar, const unsigned int version) {
        NOT_IMPLEMENTED_ABORT
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER();
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
    std::string flags ="      ";
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
    }
    o << " " << flags;
    o << " sqNorm=" << nd.squareNorm;
    if (nd.hasCoefs()) {
        o << " Coefs={";
        o << nd.getCoefs()[0] << ", " <<
                nd.getCoefs()[nd.getNCoefs() - 1] << "}";
    }
    return o;
}

#endif /* MWNODE_H_ */
