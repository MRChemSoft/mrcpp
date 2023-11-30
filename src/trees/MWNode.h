/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2021 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
 *
 * This file is part of MRCPP.
 *
 * MRCPP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRCPP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRCPP.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRCPP, see:
 * <https://mrcpp.readthedocs.io/>
 */


#pragma once

#include <Eigen/Core>

#include "MRCPP/macros.h"
#include "utils/omp_utils.h"

#include "HilbertPath.h"
#include "MWTree.h"
#include "NodeIndex.h"

namespace mrcpp {

/** @class MWNode
 *
 * @brief Base class for Multiwavelet nodes
 *
 * @details A MWNode will contain the scaling and wavelet coefficients
 * to represent functions or operators within a Multiwavelet
 * framework. The nodes are in multidimensional. The dimensionality is
 * set thoucgh the template parameter D=1,2,3. In addition to the
 * coefficients the node contains metadata such as the scale, the
 * translation index, the norm, pointers to parent node and child
 * nodes, pointer to the corresponding MWTree etc... See member and
 * data descriptions for details.
 *
 */
template <int D> class MWNode {
public:
    MWNode(const MWNode<D> &node, bool allocCoef = true, bool SetCoef = true);
    MWNode<D> &operator=(const MWNode<D> &node) = delete;
    virtual ~MWNode();

    int getKp1() const { return getMWTree().getKp1(); }
    int getKp1_d() const { return getMWTree().getKp1_d(); }
    int getOrder() const { return getMWTree().getOrder(); }
    int getScalingType() const { return getMWTree().getMRA().getScalingBasis().getScalingType(); }
    int getTDim() const { return (1 << D); }
    int getDepth() const { return getNodeIndex().getScale() - getMWTree().getRootScale(); }
    int getScale() const { return getNodeIndex().getScale(); }
    int getNChildren() const { return (isBranchNode()) ? getTDim() : 0; }
    int getSerialIx() const { return this->serialIx; }
    void setSerialIx(int Ix) { this->serialIx = Ix; }

    const NodeIndex<D> &getNodeIndex() const { return this->nodeIndex; }
    const HilbertPath<D> &getHilbertPath() const { return this->hilbertPath; }

    Coord<D> getCenter() const;
    Coord<D> getUpperBounds() const;
    Coord<D> getLowerBounds() const;

    bool hasCoord(const Coord<D> &r) const;
    bool isCompatible(const MWNode<D> &node);
    bool isAncestor(const NodeIndex<D> &idx) const;
    bool isDecendant(const NodeIndex<D> &idx) const;

    double getSquareNorm() const { return this->squareNorm; }
    double getMaxSquareNorm() const { return (maxSquareNorm > 0.0) ? maxSquareNorm : calcScaledSquareNorm(); }
    double getMaxWSquareNorm() const { return (maxWSquareNorm > 0.0) ? maxWSquareNorm : calcScaledWSquareNorm(); }

    double getScalingNorm() const;
    virtual double getWaveletNorm() const;
    double getComponentNorm(int i) const { return this->componentNorms[i]; }

    int getNCoefs() const { return this->n_coefs; }
    void getCoefs(Eigen::VectorXd &c) const;
    void printCoefs() const;

    double *getCoefs() { return this->coefs; }
    const double *getCoefs() const { return this->coefs; }

    void getPrimitiveQuadPts(Eigen::MatrixXd &pts) const;
    void getPrimitiveChildPts(Eigen::MatrixXd &pts) const;
    void getExpandedQuadPts(Eigen::MatrixXd &pts) const;
    void getExpandedChildPts(Eigen::MatrixXd &pts) const;

    MWTree<D> &getMWTree() { return static_cast<MWTree<D> &>(*this->tree); }
    MWNode<D> &getMWParent() { return static_cast<MWNode<D> &>(*this->parent); }
    MWNode<D> &getMWChild(int i) { return static_cast<MWNode<D> &>(*this->children[i]); }

    const MWTree<D> &getMWTree() const { return static_cast<const MWTree<D> &>(*this->tree); }
    const MWNode<D> &getMWParent() const { return static_cast<const MWNode<D> &>(*this->parent); }
    const MWNode<D> &getMWChild(int i) const { return static_cast<const MWNode<D> &>(*this->children[i]); }

    void zeroCoefs();
    void setCoefBlock(int block, int block_size, const double *c);
    void addCoefBlock(int block, int block_size, const double *c);
    void zeroCoefBlock(int block, int block_size);
    void attachCoefs(double *coefs);

    void calcNorms();
    void zeroNorms();
    void clearNorms();

    virtual void createChildren(bool coefs);
    virtual void genChildren();
    virtual void genParent();
    virtual void deleteChildren();
    virtual void deleteParent();

    virtual void cvTransform(int kind);
    virtual void mwTransform(int kind);

    double getNodeNorm(const NodeIndex<D> &idx) const;

    bool hasParent() const { return (parent != nullptr) ? true : false; }
    bool hasCoefs() const { return (this->status & FlagHasCoefs); }
    bool isEndNode() const { return (this->status & FlagEndNode); }
    bool isGenNode() const { return (this->status & FlagGenNode); }
    bool isRootNode() const { return (this->status & FlagRootNode); }
    bool isLeafNode() const { return not(this->status & FlagBranchNode); }
    bool isAllocated() const { return (this->status & FlagAllocated); }
    bool isBranchNode() const { return (this->status & FlagBranchNode); }
    bool isLooseNode() const { return (this->status & FlagLooseNode); }
    bool checkStatus(unsigned char mask) const { return (mask == (this->status & mask)); }

    void setHasCoefs() { SET_BITS(status, FlagHasCoefs | FlagAllocated); }
    void setIsEndNode() { SET_BITS(status, FlagEndNode); }
    void setIsGenNode() { SET_BITS(status, FlagGenNode); }
    void setIsRootNode() { SET_BITS(status, FlagRootNode); }
    void setIsLeafNode() { CLEAR_BITS(status, FlagBranchNode); }
    void setIsAllocated() { SET_BITS(status, FlagAllocated); }
    void setIsBranchNode() { SET_BITS(status, FlagBranchNode); }
    void setIsLooseNode() { SET_BITS(status, FlagLooseNode); }
    void clearHasCoefs() { CLEAR_BITS(status, FlagHasCoefs); }
    void clearIsEndNode() { CLEAR_BITS(status, FlagEndNode); }
    void clearIsGenNode() { CLEAR_BITS(status, FlagGenNode); }
    void clearIsRootNode() { CLEAR_BITS(status, FlagRootNode); }
    void clearIsAllocated() { CLEAR_BITS(status, FlagAllocated); }

    friend std::ostream &operator<<(std::ostream &o, const MWNode<D> &nd) { return nd.print(o); }

    friend class TreeBuilder<D>;
    friend class MultiplicationCalculator<D>;
    friend class NodeAllocator<D>;
    friend class MWTree<D>;
    friend class FunctionTree<D>;
    friend class OperatorTree;
    friend class FunctionNode<D>;
    friend class OperatorNode;
    friend class DerivativeCalculator<D>;

protected:
    MWTree<D> *tree{nullptr};    ///< Tree the node belongs to
    MWNode<D> *parent{nullptr};  ///< Parent node
    MWNode<D> *children[1 << D]; ///< 2^D children

    double squareNorm{-1.0};       ///< Squared norm of all 2^D (k+1)^D coefficients
    double componentNorms[1 << D]; ///< Squared norms of the separeted 2^D components
    double maxSquareNorm{-1.0};    ///< Largest squared norm among itself and descendants.
    double maxWSquareNorm{-1.0};   ///< Largest wavelet squared norm among itself and descendants.
                                   ///< NB: must be set before used.
    double *coefs{nullptr};     ///< the 2^D (k+1)^D MW coefficients
                                ///< For example, in case of a one dimensional function \f$ f \f$
                                ///< this array equals \f$ s_0, \ldots, s_k, d_0, \ldots, d_k \f$,
                                ///< where scaling coefficients \f$ s_j = s_{jl}^n(f) \f$
                                ///< and wavelet coefficients \f$ d_j = d_{jl}^n(f) \f$.
                                ///< Here \f$ n, l \f$ are unique for every node.
    int n_coefs{0};

    int serialIx{-1};       ///< index in serial Tree
    int parentSerialIx{-1}; ///< index of parent in serial Tree, or -1 for roots
    int childSerialIx{-1};  ///< index of first child in serial Tree, or -1 for leafnodes/endnodes

    NodeIndex<D> nodeIndex;     ///< Scale and translation of the node
    HilbertPath<D> hilbertPath; ///< To be documented

    MWNode();
    MWNode(MWTree<D> *tree, int rIdx);
    MWNode(MWTree<D> *tree, const NodeIndex<D> &idx);
    MWNode(MWNode<D> *parent, int cIdx);
    virtual void dealloc();

    bool crop(double prec, double splitFac, bool absPrec);

    void initNodeLock() { MRCPP_INIT_OMP_LOCK(); }
    virtual void allocCoefs(int n_blocks, int block_size);
    virtual void freeCoefs();

    void setMaxSquareNorm();
    void resetMaxSquareNorm();
    double calcScaledSquareNorm() const { return std::pow(2.0, D * getScale()) * getSquareNorm(); }
    double calcScaledWSquareNorm() const { return std::pow(2.0, D * getScale()) * getWaveletNorm(); }
    virtual double calcComponentNorm(int i) const;

    virtual void reCompress();
    virtual void giveChildrenCoefs(bool overwrite = true);
    virtual void giveChildCoefs(int cIdx, bool overwrite = true);
    virtual void giveParentCoefs(bool overwrite = true);
    virtual void copyCoefsFromChildren();

    int getChildIndex(const NodeIndex<D> &nIdx) const;
    int getChildIndex(const Coord<D> &r) const;

    bool diffBranch(const MWNode<D> &rhs) const;

    MWNode<D> *retrieveNode(const Coord<D> &r, int depth);
    MWNode<D> *retrieveNode(const NodeIndex<D> &idx);
    MWNode<D> *retrieveParent(const NodeIndex<D> &idx);

    const MWNode<D> *retrieveNodeNoGen(const NodeIndex<D> &idx) const;
    MWNode<D> *retrieveNodeNoGen(const NodeIndex<D> &idx);

    const MWNode<D> *retrieveNodeOrEndNode(const Coord<D> &r, int depth) const;
    MWNode<D> *retrieveNodeOrEndNode(const Coord<D> &r, int depth);

    const MWNode<D> *retrieveNodeOrEndNode(const NodeIndex<D> &idx) const;
    MWNode<D> *retrieveNodeOrEndNode(const NodeIndex<D> &idx);

    void threadSafeGenChildren();
    void deleteGenerated();

    virtual std::ostream &print(std::ostream &o) const;

    static const unsigned char FlagBranchNode = B8(00000001);
    static const unsigned char FlagGenNode = B8(00000010);
    static const unsigned char FlagHasCoefs = B8(00000100);
    static const unsigned char FlagAllocated = B8(00001000);
    static const unsigned char FlagEndNode = B8(00010000);
    static const unsigned char FlagRootNode = B8(00100000);
    static const unsigned char FlagLooseNode = B8(01000000);

private:
    unsigned char status{0};

#ifdef MRCPP_HAS_OMP
    omp_lock_t omp_lock;
#endif
};

} // namespace mrcpp
