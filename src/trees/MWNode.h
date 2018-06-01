/**
 *  Simple n-dimensional node
 */

#pragma once

#pragma GCC system_header
#include <Eigen/Core>

#include "utils/omp_utils.h"
#include "macros.h"

#include "MWTree.h"
#include "NodeIndex.h"
#include "HilbertPath.h"

namespace mrcpp {

template<int D>
class MWNode {
public:
    MWNode(const MWNode<D> &node);
    virtual ~MWNode();

    int getKp1() const { return getMWTree().getKp1(); }
    int getKp1_d() const { return getMWTree().getKp1_d(); }
    int getOrder() const { return getMWTree().getOrder(); }
    int getScalingType() const { return getMWTree().getMRA().getScalingBasis().getScalingType(); }
    int getTDim() const { return getMWTree().getTDim(); }
    int getDepth() const { return getNodeIndex().getScale()-getMWTree().getRootScale(); }
    int getScale() const { return getNodeIndex().getScale(); }
    int getNChildren() const { if (isBranchNode()) return getTDim(); return 0; }
    int getSerialIx() const { return this->serialIx; }
    void setSerialIx(int Ix) { this->serialIx = Ix; }
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

    void getPrimitiveQuadPts(Eigen::MatrixXd &pts) const;
    void getPrimitiveChildPts(Eigen::MatrixXd &pts) const;
    void getExpandedQuadPts(Eigen::MatrixXd &pts) const;
    void getExpandedChildPts(Eigen::MatrixXd &pts) const;

    MWTree<D>& getMWTree() { return static_cast<MWTree<D> &>(*this->tree); }
    MWNode<D>& getMWParent() { return static_cast<MWNode<D> &>(*this->parent); }
    MWNode<D>& getMWChild(int i) { return static_cast<MWNode<D> &>(*this->children[i]); }

    const MWTree<D>& getMWTree() const { return static_cast<const MWTree<D> &>(*this->tree); }
    const MWNode<D>& getMWParent() const { return static_cast<const MWNode<D> &>(*this->parent); }
    const MWNode<D>& getMWChild(int i) const { return static_cast<const MWNode<D> &>(*this->children[i]); }

    void zeroCoefs();
    void setCoefBlock(int block, int block_size, const double *c);
    void addCoefBlock(int block, int block_size, const double *c);
    void zeroCoefBlock(int block, int block_size);

    void calcNorms();
    void zeroNorms();
    void clearNorms();

    virtual void createChildren();
    virtual void genChildren();
    virtual void deleteChildren();

    virtual void cvTransform(int kind);
    virtual void mwTransform(int kind);

    bool splitCheck(double prec, double splitFac, bool absPrec) const;

    void setHasCoefs() { SET_BITS(status, FlagHasCoefs | FlagAllocated); }
    void setIsEndNode() { SET_BITS(status, FlagEndNode); }
    void setIsGenNode() { SET_BITS(status, FlagGenNode); }
    void setIsRootNode() { SET_BITS(status, FlagRootNode); }
    void setIsLeafNode() { CLEAR_BITS(status, FlagBranchNode); }
    void setIsAllocated() { SET_BITS(status, FlagAllocated); }
    void setIsBranchNode() { SET_BITS(status, FlagBranchNode); }
    void setIsLooseNode() { SET_BITS(status, FlagLooseNode); }
    void clearHasCoefs() { CLEAR_BITS(status, FlagHasCoefs);}
    void clearIsEndNode() { CLEAR_BITS(status, FlagEndNode); }
    void clearIsGenNode() { CLEAR_BITS(status, FlagGenNode); }
    void clearIsRootNode() { CLEAR_BITS(status, FlagRootNode); }
    void clearIsAllocated() { CLEAR_BITS(status, FlagAllocated); }

    friend std::ostream& operator<<(std::ostream &o, const MWNode<D> &nd) { return nd.print(o); }

    friend class TreeBuilder<D>;
    friend class MultiplicationCalculator<D>;
    friend class SerialFunctionTree<D>;
    friend class SerialOperatorTree;
    friend class MWTree<D>;
    friend class FunctionTree<D>;
    friend class OperatorTree;

protected:
    MWTree<D> *tree;
    MWNode<D> *parent;	    ///< Parent node
    MWNode<D> *children[1<<D];    ///< 2^D children

    double squareNorm;
    double componentNorms[1<<D]; ///< 2^D components

    double *coefs;
    int n_coefs;

    int lockX;          //manual lock tag (avoiding omp set/unset)
    int serialIx;       //index in serial Tree
    int parentSerialIx; //index of parent in serial Tree, or -1 for roots
    int childSerialIx;  //index of first child in serial Tree, or -1 for leafnodes/endnodes

    NodeIndex<D> nodeIndex;
    HilbertPath<D> hilbertPath;

    MWNode();
    virtual void dealloc();

    bool crop(double prec, double splitFac, bool absPrec);
    double getScaleFactor(double splitFac, bool absPrec) const;

    virtual void allocCoefs(int n_blocks, int block_size);
    virtual void freeCoefs();

    virtual double calcComponentNorm(int i) const;

    bool crop(double prec, NodeIndexSet *cropIdx = 0);

    virtual void reCompress();
    virtual void giveChildrenCoefs(bool overwrite = true);
    virtual void copyCoefsFromChildren();

    int getChildIndex(const NodeIndex<D> &nIdx) const;
    int getChildIndex(const double *r) const;

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

    void threadSafeGenChildren();
    void deleteGenerated();

    virtual std::ostream& print(std::ostream &o) const;

    static const unsigned char FlagBranchNode = B8(00000001);
    static const unsigned char FlagGenNode    = B8(00000010);
    static const unsigned char FlagHasCoefs   = B8(00000100);
    static const unsigned char FlagAllocated  = B8(00001000);
    static const unsigned char FlagEndNode    = B8(00010000);
    static const unsigned char FlagRootNode   = B8(00100000);
    static const unsigned char FlagLooseNode  = B8(01000000);

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

}
