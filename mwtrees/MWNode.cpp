/**
 *  Simple n-dimensional node
 *
 *  Created on: Oct 12, 2009
 *      Author: jonas
 */

#include "MWNode.h"
#include "MWTree.h"
#include "ProjectedNode.h"
#include "MathUtils.h"
#include "QuadratureCache.h"
#include "GenNode.h"

using namespace std;
using namespace Eigen;

/** MWNode rootnode constructor.
  * Creates an empty rootnode given its tree and node index */
template<int D>
MWNode<D>::MWNode(MWTree<D> &t, const NodeIndex<D> &nIdx)
        : tree(&t),
          parent(0),
          nodeIndex(nIdx),
          hilbertPath(),
          squareNorm(-1.0),
          status(0),
          coefs(0) {
    clearNorms();
    this->coefs = &this->coefvec;
    this->tree->incrementNodeCount(getScale());
    for (int i = 0; i < getTDim(); i++) {
        this->children[i] = 0;
    }
    setIsLeafNode();
    setIsRootNode();

#ifdef OPENMP
    omp_init_lock(&node_lock);
#endif
}

/** MWNode constructor.
  * Creates an empty node given its parent and child index */
template<int D>
MWNode<D>::MWNode(MWNode<D> &p, int cIdx)
        : tree(p.tree),
          parent(&p),
          nodeIndex(p.getNodeIndex(), cIdx),
          hilbertPath(p.getHilbertPath(), cIdx),
          squareNorm(-1.0),
          status(0),
          coefs(0) {
    clearNorms();
    this->coefs = &this->coefvec;
    this->tree->incrementNodeCount(getScale());
    for (int i = 0; i < getTDim(); i++) {
        this->children[i] = 0;
    }
    setIsLeafNode();

#ifdef OPENMP
    omp_init_lock(&node_lock);
#endif
}

/** Detatched node */
template<int D>
MWNode<D>::MWNode(const MWNode<D> &n)
        : tree(n.tree),
          parent(0),
          nodeIndex(n.getNodeIndex()),
          hilbertPath(n.getHilbertPath()),
          squareNorm(-1.0),
          status(0),
          coefs(0) {

    this->coefs = &this->coefvec;
    allocCoefs(this->getTDim());
    if (n.hasCoefs()) {
        setCoefs(n.getCoefs());
    } else {
        zeroCoefs();
    }
    for (int i = 0; i < getTDim(); i++) {
        this->children[i] = 0;
    }
    setIsLeafNode();
    setIsLooseNode();

#ifdef OPENMP
    omp_init_lock(&node_lock);
#endif
}

/** MWNode destructor.
  * Recursive deallocation of a node and all its decendants */
template<int D>
MWNode<D>::~MWNode() {
    if (this->isBranchNode()) {
        deleteChildren();
    }
    if (not this->isLooseNode()) {
        this->tree->decrementNodeCount(getScale());
	if (this->NodeRank<0){println(0, "Node has no rank! " << this->NodeRank);
	}else{        this->tree->allocator->DeAllocNodes(this->NodeRank);}
    } else {
        freeCoefs();
    }
    *((double**) ((void*)&(this->coefvec))) = NULL;
    *((int*) (((void*)&(this->coefvec))+8)) = 0;
#ifdef OPENMP
    omp_destroy_lock(&node_lock);
#endif
}

/** Allocate the coefs vector. */
template<int D>
void MWNode<D>::allocCoefs(int nBlocks) {
    if (this->isAllocated()) MSG_FATAL("Coefs already allocated");

    int nCoefs = nBlocks * this->getKp1_d();
    if(*((int*) (((void*)&(this->coefvec))+8)) != 0)println(0, this->GenNodeCoeffIx<<" Node already allocated!   "<<this->NodeCoeffIx);
    if (ProjectedNode<D> *node = dynamic_cast<ProjectedNode<D> *>(this)) {
      //      this->coefs = this->tree->allocator->allocCoeff(nBlocks);
      *((double**) ((void*)&(this->coefvec))) = this->tree->allocator->allocCoeff(nBlocks);
      *((int*) (((void*)&(this->coefvec))+8)) = this->tree->allocator->sizeNodeCoeff/sizeof(double);
      this->NodeCoeffIx = this->tree->allocator->nNodesCoeff;
    } else if (GenNode<D> *node = dynamic_cast<GenNode<D> *>(this)) {
	//      this->coefs = this->tree->allocator->allocGenCoeff(nBlocks);
      //*((double**) ((void*)&(this->coefvec))) = this->tree->allocator->allocGenCoeff(nBlocks);
      //*((int*) (((void*)&(this->coefvec))+8)) = this->tree->allocator->sizeGenNodeCoeff/sizeof(double);
       // this->GenNodeCoeffIx = this->tree->allocator->nGenNodesCoeff;
      *((double**) ((void*)&(this->coefvec))) = this->tree->allocator->allocCoeff(nBlocks);
      *((int*) (((void*)&(this->coefvec))+8)) = this->tree->allocator->sizeNodeCoeff/sizeof(double);
      this->GenNodeCoeffIx = this->tree->allocator->nNodesCoeff;
      this->NodeCoeffIx = this->tree->allocator->nNodesCoeff;
   } else{//default is to allocate on ProjectedNode stack
	//this->coefs = this->tree->allocator->allocCoeff(nBlocks);
      *((double**) ((void*)&(this->coefvec))) = this->tree->allocator->allocCoeff(nBlocks);
      *((int*) (((void*)&(this->coefvec))+8)) = this->tree->allocator->sizeNodeCoeff/sizeof(double);
      this->NodeCoeffIx = this->tree->allocator->nNodesCoeff;
    }

    this->setIsAllocated();
    this->clearHasCoefs();
}

/** Deallocation of coefficients. */
template<int D>
void MWNode<D>::freeCoefs() {
    if (not this->isAllocated()) MSG_FATAL("Coefs not allocated");
    //delete this->coefs;
    //this->coefs->~VectorXd();
    if (ProjectedNode<D> *node = dynamic_cast<ProjectedNode<D> *>(this)) {
      this->tree->allocator->DeAllocCoeff(this->NodeCoeffIx);
    } else if (GenNode<D> *node = dynamic_cast<GenNode<D> *>(this)) {
      //this->tree->allocator->DeAllocGenCoeff(this->GenNodeCoeffIx);
      this->tree->allocator->DeAllocCoeff(this->NodeCoeffIx);
    } else{//default is to (de)allocate on ProjectedNode stack
      this->tree->allocator->DeAllocCoeff(this->NodeCoeffIx);
    }

    
    this->coefs = 0;
    this->clearHasCoefs();
    this->clearIsAllocated();
}

template<int D>
void MWNode<D>::zeroCoefs() {
    if (not this->isAllocated()) MSG_FATAL("Coefs not allocated");
    this->coefs->setZero();
    this->zeroNorms();
    this->setHasCoefs();
}

/** Set coefficients of node.
  *
  * Copies the argument vector to the coefficient vector of the node. */
template<int D>
void MWNode<D>::setCoefs(const Eigen::VectorXd &c) {
    if (not this->isAllocated()) MSG_FATAL("Coefs not allocated");

    int nNew = c.size();
    int nEmpty = this->getNCoefs() - nNew;
    if (nEmpty < 0) {
        MSG_FATAL("Size mismatch");
    } else {
        this->coefs->segment(nNew, nEmpty).setZero();
    }
    this->coefs->segment(0, nNew) = c;
    this->setHasCoefs();
    this->calcNorms();
}

template<int D>
void MWNode<D>::giveChildrenCoefs(bool overwrite) {
    assert(this->isBranchNode());
    if (not this->hasCoefs()) MSG_FATAL("No coefficients!");

    MWNode<D> copy(*this);
    copy.mwTransform(Reconstruction);
    const VectorXd &c = copy.getCoefs();

    int kp1_d = this->getKp1_d();
    for (int i = 0; i < this->getTDim(); i++) {
        MWNode<D> &child = this->getMWChild(i);
        if (overwrite) {
            child.setCoefs(c.segment(i*kp1_d, kp1_d));
        } else if (child.hasCoefs()) {
            child.getCoefs().segment(0, kp1_d) += c.segment(i*kp1_d, kp1_d);
        } else {
            MSG_FATAL("Child has no coefs");
        }
        child.calcNorms();
    }
}

/** Takes the scaling coefficients of the children and stores them consecutively
  * in the  given vector. */
template<int D>
void MWNode<D>::copyCoefsFromChildren(VectorXd &c) {
    int kp1_d = this->getKp1_d();
    for (int i = 0; i < this->getTDim(); i++) {
        MWNode<D> &child = getMWChild(i);
        if (child.hasCoefs()) {
            VectorXd &cc = child.getCoefs();
            c.segment(i * kp1_d, kp1_d) = cc.segment(0, kp1_d);
	   //for (int j = 0; j < kp1_d; j++) c[j+i * kp1_d]=child.getCoefs()[j];
        } else {
	    c.segment(i * kp1_d, kp1_d).setZero();
	  //   for (int j = 0; j < kp1_d; j++) c[j+i * kp1_d]=0.0;
        }
    }
}


/** Coefficient-Value transform
  *
  * This routine transforms the scaling coefficients of the node to the
  * function values in the corresponding quadrature roots (of its children).
  * Input parameter = forward: coef->value.
  * Input parameter = backward: value->coef.
  *
  * NOTE: this routine assumes a 0/1 (scaling on children 0 and 1)
  *       representation, in oppose to s/d (scaling and wavelet). */
template<int D>
void MWNode<D>::cvTransform(int operation) {
    const ScalingBasis &sf = this->getMWTree().getMRA().getScalingBasis();
    if (sf.getScalingType() != Interpol) {
        NOT_IMPLEMENTED_ABORT;
    }

    int quadratureOrder = sf.getQuadratureOrder();
    getQuadratureCache(qc);

    double two_scale = pow(2.0, this->getScale() + 1);
    VectorXd modWeights = qc.getWeights(quadratureOrder);
    if (operation == Forward) {
        modWeights = modWeights.array().inverse();
        modWeights *= two_scale;
        modWeights = modWeights.array().sqrt();
    } else if (operation == Backward) {
        modWeights *= 1.0/two_scale;
        modWeights = modWeights.array().sqrt();
    } else {
        MSG_FATAL("Invalid operation");
    }

    VectorXd &coefs = this->getCoefs();

    int kp1 = this->getKp1();
    int kp1_d = this->getKp1_d();
    int kp1_p[D];
    for (int d = 0; d < D; d++) {
        kp1_p[d] = MathUtils::ipow(kp1, d);
    }

    for (int m = 0; m < this->getTDim(); m++) {
        for (int p = 0; p < D; p++) {
            int n = 0;
            for (int i = 0; i < kp1_p[D - p - 1]; i++) {
                for (int j = 0; j < kp1; j++) {
                    for (int k = 0; k < kp1_p[p]; k++) {
                        coefs[m * kp1_d + n] *= modWeights[j];
                        n++;
                    }
                }
            }
        }
    }
}

/** Multiwavelet transform: fast version
  *
  * Application of the filters on one node to pass from a 0/1 (scaling
  * on children 0 and 1) representation to an s/d (scaling and
  * wavelet) representation. Bit manipulation is used in order to
  * determine the correct filters and whether to apply them or just
  * pass to the next couple of indexes. The starting coefficients are
  * preserved until the application is terminated, then they are
  * overwritten. With minor modifications this code can also be used
  * for the inverse mw transform (just use the transpose filters) or
  * for the application of an operator (using A, B, C and T parts of an
  * operator instead of G1, G0, H1, H0). This is the version where the
  * three directions are operated one after the other. Although this
  * is formally faster than the other algorithm, the separation of the
  * three dimensions prevent the possibility to use the norm of the
  * operator in order to discard a priori negligible contributions.

  * Luca Frediani, August 2006
  * C++ version: Jonas Juselius, September 2009 */
template<int D>
void MWNode<D>::orig_mwTransform(int operation) {
    int kp1 = this->getKp1();
    int kp1_dm1 = MathUtils::ipow(kp1, D - 1);
    int kp1_d = this->getKp1_d();
    const MWFilter &filter = getMWTree().getMRA().getFilter();
    VectorXd &result = getMWTree().getTmpMWCoefs();
    double overwrite = 0.0;

    for (int i = 0; i < D; i++) {
        int mask = 1 << i;
        for (int gt = 0; gt < this->getTDim(); gt++) {
            double *out = result.data() + gt * kp1_d;
            for (int ft = 0; ft < this->getTDim(); ft++) {
                /* Operate in direction i only if the bits along other
                 * directions are identical. The bit of the direction we
                 * operate on determines the appropriate filter/operator */
                if ((gt | mask) == (ft | mask)) {
                    double *in = this->coefs->data() + ft * kp1_d;
                    int fIdx = 2 * ((gt >> i) & 1) + ((ft >> i) & 1);
                    const MatrixXd &oper = filter.getSubFilter(fIdx, operation);
                    MathUtils::applyFilter(out, in, oper, kp1, kp1_dm1, overwrite);
                    overwrite = 1.0;
                }
            }
            overwrite = 0.0;
        }
        this->coefs->swap(result);
    }
}

/** Multiwavelet transform: fast version
  *
  * Application of the filters on one node to pass from a 0/1 (scaling
  * on children 0 and 1) representation to an s/d (scaling and
  * wavelet) representation. Bit manipulation is used in order to
  * determine the correct filters and whether to apply them or just
  * pass to the next couple of indexes. The starting coefficients are
  * preserved until the application is terminated, then they are
  * overwritten. With minor modifications this code can also be used
  * for the inverse mw transform (just use the transpose filters) or
  * for the application of an operator (using A, B, C and T parts of an
  * operator instead of G1, G0, H1, H0). This is the version where the
  * three directions are operated one after the other. Although this
  * is formally faster than the other algorithm, the separation of the
  * three dimensions prevent the possibility to use the norm of the
  * operator in order to discard a priori negligible contributions.

  * Luca Frediani, August 2006
  * C++ version: Jonas Juselius, September 2009 */
template<int D>
void MWNode<D>::mwTransform(int operation) {
    int kp1 = this->getKp1();
    int kp1_dm1 = MathUtils::ipow(kp1, D - 1);
    int kp1_d = this->getKp1_d();
    const MWFilter &filter = getMWTree().getMRA().getFilter();
    VectorXd &result = getMWTree().getTmpMWCoefs();
    double overwrite = 0.0;
    double *data_p1,*data_p2,*data_p3;

    data_p1 = this->coefs->data();
    data_p2 = result.data();
        
    for (int i = 0; i < D; i++) {
      int mask = 1 << i;
      for (int gt = 0; gt < this->getTDim(); gt++) {
	double *out =  data_p2 + gt * kp1_d;
	for (int ft = 0; ft < this->getTDim(); ft++) {
	  /* Operate in direction i only if the bits along other
	   * directions are identical. The bit of the direction we
	   * operate on determines the appropriate filter/operator */
	  if ((gt | mask) == (ft | mask)) {
	    double *in = data_p1 + ft * kp1_d;
	    int fIdx = 2 * ((gt >> i) & 1) + ((ft >> i) & 1);
	    const MatrixXd &oper = filter.getSubFilter(fIdx, operation);
	    MathUtils::applyFilter(out, in, oper, kp1, kp1_dm1, overwrite);
	    overwrite = 1.0;
	  }
	}
	overwrite = 0.0;
      }
      //swap pointers (but data is not copied!)
      data_p3=data_p1;
      data_p1=data_p2;
      data_p2=data_p3;
    }
        
    //if the loop is traversed an odd number of time the data ends in the wrong buffer
    //(but now data is in cache, so this is very fast)
    if(D%2)for(int i=0; i<kp1_d*this->getTDim(); i++) this->coefvec[i] = data_p1[i];
    
}

/** Set all norms to Undefined. */
template<int D>
void MWNode<D>::clearNorms() {
    this->squareNorm = -1.0;
    for (int i = 0; i < this->getTDim(); i++) {
        this->componentNorms[i] = -1.0;
    }
}

/** Set all norms to zero. */
template<int D>
void MWNode<D>::zeroNorms() {
    this->squareNorm = 0.0;
    for (int i = 0; i < this->getTDim(); i++) {
        this->componentNorms[i] = 0.0;
    }
}

/** Calculate and store square norm and component norms, if allocated. */
template<int D>
void MWNode<D>::calcNorms() {
    this->squareNorm = 0.0;
    for (int i = 0; i < this->getTDim(); i++) {
        double norm_i = calcComponentNorm(i);
        this->componentNorms[i] = norm_i;
        this->squareNorm += norm_i*norm_i;
    }
}

/** Calculate and return the squared scaling norm. */
template<int D>
double MWNode<D>::getScalingNorm() const {
    double sNorm = this->getComponentNorm(0);
    if (sNorm >= 0.0) {
        return sNorm*sNorm;
    } else {
        return -1.0;
    }
}

/** Calculate and return the squared wavelet norm. */
template<int D>
double MWNode<D>::getWaveletNorm() const {
    double wNorm = 0.0;
    for (int i = 1; i < this->getTDim(); i++) {
        double norm_i = this->getComponentNorm(i);
        if (norm_i >= 0.0) {
            wNorm += norm_i*norm_i;
        } else {
            wNorm = -1.0;
        }
    }
    return wNorm;
}

/** Calculate the norm of one component (NOT the squared norm!). */
template<int D>
double MWNode<D>::calcComponentNorm(int i) const {
    assert(this->isAllocated());
    assert(this->hasCoefs());
    const VectorXd &c = this->getCoefs();
    int kp1_d = this->getKp1_d();
    return c.segment(i*kp1_d, kp1_d).norm();
}

template<int D>
double MWNode<D>::estimateError(bool absPrec) {
    NOT_IMPLEMENTED_ABORT;
    //    if (this->isForeign()) {
    //        return 0.0;
    //    }
    //    if (this->isCommon() and this->tree->getRankId() != 0) {
    //        return 0.0;
    //    }
    //    double tNorm = 1.0;
    //    if (not absPrec) {
    //        tNorm = sqrt(getMWTree().getSquareNorm());
    //    }

    //    int n = this->getScale();
    //    double expo = (1.0 * (n + 1));
    //    double scaleFactor = max(2.0* MachinePrec, pow(2.0, -expo));
    //    double wNorm = this->calcWaveletNorm();
    //    double error = scaleFactor * wNorm / tNorm;
    //    return error*error;
}

/** Update the coefficients of the node by a mw transform of the scaling
  * coefficients of the children. Option to overwrite or add up existing
  * coefficients. */
template<int D>
void MWNode<D>::reCompress(bool overwrite) {
    if ((not this->isGenNode()) and this->isBranchNode()) {
        if (not this->isAllocated()) MSG_FATAL("Coefs not allocated");
        if (overwrite) {
            copyCoefsFromChildren(*this->coefs);
	    mwTransform(Compression);
	    //            int Children_Stride = this->getMWTree().allocator->sizeGenNodeCoeff/sizeof(double);
            //            this->getMWTree().allocator->S_mwTransformBack(&(this->children[0]->getCoefs()[0]),&(this->getCoefs()[0]),Children_Stride);

        } else {
            // Check optimization
            NOT_IMPLEMENTED_ABORT;
            MatrixXd tmp = getCoefs();
            copyCoefsFromChildren(*this->coefs);
            mwTransform(Compression);
            getCoefs() += tmp;
        }
        this->setHasCoefs();
        calcNorms();
    }
}

/** Recurse down until an EndNode is found, and then crop children with
  * too high precision. */
template<int D>
bool MWNode<D>::crop(double prec, NodeIndexSet *cropIdx) {
    NOT_IMPLEMENTED_ABORT;
    //    if (this->isEndNode()) {
    //        return true;
    //    } else {
    //        for (int i = 0; i < this->tDim; i++) {
    //            MWNode<D> &child = *this->children[i];
    //            if (child.cropChildren(prec, cropIdx)) {
    //                if (not this->isForeign()) {
    //                    if (this->splitCheck(prec) == false) {
    //                        if (cropIdx != 0) {
    //                            cropIdx->insert(&this->getNodeIndex());
    //                        } else {
    //                            this->deleteChildren();
    //                        }
    //                        return true;
    //                    }
    //                }
    //            }
    //        }
    //    }
    //    return false;
}

template<int D>
void MWNode<D>::createChildren() {
    if (this->isBranchNode()) MSG_FATAL("Node already has children");
    for (int cIdx = 0; cIdx < getTDim(); cIdx++) {
        createChild(cIdx);
    }
    this->setIsBranchNode();
}

/** Recursive deallocation of children and all their decendants.
  * Leaves node as LeafNode and children[] as null pointer. */
template<int D> void MWNode<D>::deleteChildren() {
     for (int cIdx = 0; cIdx < getTDim(); cIdx++) {
         if (this->children[cIdx] != 0) {
             // ProjectedNode<D> *node = static_cast<ProjectedNode<D> *>(this->children[cIdx]);
             // delete node;             // node->~ProjectedNode();
             if (ProjectedNode<D> *node = dynamic_cast<ProjectedNode<D> *>(this->children[cIdx])) {
                 node->~ProjectedNode();
             } else if (GenNode<D> *node = dynamic_cast<GenNode<D> *>(this->children[cIdx])) {
	       node->~GenNode();
             }
	     //             delete this->children[cIdx];
             this->children[cIdx] = 0;
         }
     } 
    this->setIsLeafNode();
 }

template<int D>
void MWNode<D>::genChildren() {
    if (this->isBranchNode()) MSG_FATAL("Node already has children");
    int nChildren = this->getTDim();
    for (int cIdx = 0; cIdx < nChildren; cIdx++) {
        genChild(cIdx);
    }
    this->setIsBranchNode();
}

/** Clear coefficients of generated nodes.
  *
  * The node structure is kept, only the coefficients are cleared. */
template<int D>
void MWNode<D>::clearGenerated() {
    if (this->isBranchNode()) {
        assert(this->children != 0);
        for (int cIdx = 0; cIdx < this->getTDim(); cIdx++) {
            if (this->children[cIdx] != 0) {
                this->getMWChild(cIdx).clearGenerated();
            }
        }
    }
}

template<int D>
void MWNode<D>::deleteGenerated() {
    if (this->isBranchNode()) {
        if (this->isEndNode()) {
            this->deleteChildren();
        } else {
            for (int cIdx = 0; cIdx < getTDim(); cIdx++) {
                assert(this->children[cIdx] != 0);
                this->getMWChild(cIdx).deleteGenerated();
            }
        }
    }
}

template<int D>
void MWNode<D>::getCenter(double *r) const {
    NOT_IMPLEMENTED_ABORT;
    //    assert(r != 0);
    //    double sFac = pow(2.0, -getScale());
    //    for (int d = 0; d < D; d++) {
    //        double l = (double) getTranslation()[d];
    //        r[d] = sFac*(l + 0.5);
    //    }
}

template<int D>
void MWNode<D>::getBounds(double *lb, double *ub) const {
    int n = getScale();
    double p = pow(2.0, -n);
    const int *l = getTranslation();
    for (int i = 0; i < D; i++) {
        lb[i] = p * l[i];
        ub[i] = p * (l[i] + 1);
    }
}

/** Routine to find the path along the tree.
  *
  * Given the translation indices at the final scale, computes the child m
  * to be followed at the current scale in oder to get to the requested
  * node at the final scale. The result is the index of the child needed.
  * The index is obtained by bit manipulation of of the translation indices. */
template<int D>
int MWNode<D>::getChildIndex(const NodeIndex<D> &nIdx) const {
    assert(isAncestor(nIdx));
    int cIdx = 0;
    int diffScale = nIdx.getScale() - getScale() - 1;
    assert(diffScale >= 0);
    for (int d = 0; d < D; d++) {
        int bit = (nIdx.getTranslation()[d] >> (diffScale)) & 1;
        cIdx = cIdx + (bit << d);
    }
    assert(cIdx >= 0);
    assert(cIdx < getTDim());
    return cIdx;
}

/** Routine to find the path along the tree.
  *
  * Given a point in space, determines which child should be followed
  * to get to the corresponding terminal node. */
template<int D>
int MWNode<D>::getChildIndex(const double *r) const {
    assert(hasCoord(r));
    int cIdx = 0;
    double sFac = pow(2.0, -getScale());
    const int *l = getTranslation();
    for (int d = 0; d < D; d++) {
        if (r[d] > sFac*(l[d] + 0.5)) {
            cIdx = cIdx + (1 << d);
        }
    }
    assert(cIdx >= 0);
    assert(cIdx < getTDim());
    return cIdx;
}

/** Const version of node retriever that NEVER generates.
  *
  * Recursive routine to find and return the node with a given NodeIndex.
  * This routine returns the appropriate ProjectedNode, or a NULL pointer if
  * the node does not exist, or if it is a GenNode. Recursion starts at at this
  * node and ASSUMES the requested node is in fact decending from this node. */
template<int D>
const MWNode<D> *MWNode<D>::retrieveNodeNoGen(const NodeIndex<D> &idx) const {
    if (getScale() == idx.getScale()) { // we're done
        assert(getNodeIndex() == idx);
        return this;
    }
    assert(this->isAncestor(idx));
    if (this->isEndNode()) { // don't return GenNodes
        return 0;
    }
    int cIdx = getChildIndex(idx);
    assert(this->children[cIdx] != 0);
    return this->children[cIdx]->retrieveNodeNoGen(idx);
}

/** Node retriever that NEVER generates.
  *
  * Recursive routine to find and return the node with a given NodeIndex.
  * This routine returns the appropriate ProjectedNode, or a NULL pointer if
  * the node does not exist, or if it is a GenNode. Recursion starts at at this
  * node and ASSUMES the requested node is in fact decending from this node. */
template<int D>
MWNode<D> *MWNode<D>::retrieveNodeNoGen(const NodeIndex<D> &idx) {
    if (getScale() == idx.getScale()) { // we're done
        assert(getNodeIndex() == idx);
        return this;
    }
    assert(this->isAncestor(idx));
    if (this->isEndNode()) { // don't return GenNodes
        return 0;
    }
    int cIdx = getChildIndex(idx);
    assert(this->children[cIdx] != 0);
    return this->children[cIdx]->retrieveNodeNoGen(idx);
}

template<int D>
const MWNode<D> *MWNode<D>::retrieveNodeOrEndNode(const double *r, int depth) const {
    if (getDepth() == depth or this->isEndNode()) {
        return this;
    }
    int cIdx = getChildIndex(r);
    assert(this->children[cIdx] != 0);
    return this->children[cIdx]->retrieveNodeOrEndNode(r, depth);
}

/** Node retriever that return requested ProjectedNode or EndNode.
  *
  * Recursive routine to find and return the node with a given NodeIndex.
  * This routine returns the appropriate ProjectedNode, or the EndNode on the
  * path to the requested node, and will never create or return GenNodes.
  * Recursion starts at at this node and ASSUMES the requested node is in fact
  * decending from this node. */
template<int D>
MWNode<D> *MWNode<D>::retrieveNodeOrEndNode(const double *r, int depth) {
    if (getDepth() == depth or this->isEndNode()) {
        return this;
    }
    int cIdx = getChildIndex(r);
    assert(this->children[cIdx] != 0);
    return this->children[cIdx]->retrieveNodeOrEndNode(r, depth);
}

template<int D>
const MWNode<D> *MWNode<D>::retrieveNodeOrEndNode(const NodeIndex<D> &idx) const {
    if (getScale() == idx.getScale()) { // we're done
        assert(getNodeIndex() == idx);
        return this;
    }
    assert(isAncestor(idx));
    // We should in principle lock before read, but it makes things slower,
    // and the EndNode status does not change (normally ;)
    if (isEndNode()) {
        return this;
    }
    int cIdx = getChildIndex(idx);
    assert(children[cIdx] != 0);
    return this->children[cIdx]->retrieveNodeOrEndNode(idx);
}

template<int D>
MWNode<D> *MWNode<D>::retrieveNodeOrEndNode(const NodeIndex<D> &idx) {
    if (getScale() == idx.getScale()) { // we're done
        assert(getNodeIndex() == idx);
        return this;
    }
    assert(isAncestor(idx));
    // We should in principle lock before read, but it makes things slower,
    // and the EndNode status does not change (normally ;)
    if (isEndNode()) {
        return this;
    }
    int cIdx = getChildIndex(idx);
    assert(children[cIdx] != 0);
    return this->children[cIdx]->retrieveNodeOrEndNode(idx);
}

/** Node retriever that ALWAYS returns the requested node.
  *
  * Recursive routine to find and return the node with a given NodeIndex.
  * This routine always returns the appropriate node, and will generate nodes
  * that does not exist. Recursion starts at at this node and ASSUMES the
  * requested node is in fact decending from this node. */
template<int D>
MWNode<D> *MWNode<D>::retrieveNode(const double *r, int depth) {
    if (depth < 0) MSG_FATAL("Invalid argument");

    if (getDepth() == depth) {
        return this;
    }
    assert(hasCoord(r));
    // If we have reached an endNode, lock if necessary, and start generating
    // NB! retrieveNode() for GenNodes behave a bit differently.
    lockNode();
    if (this->isLeafNode()) {
        genChildren();
        giveChildrenCoefs();
    }
    unlockNode();
    int cIdx = getChildIndex(r);
    assert(this->children[cIdx] != 0);
    return this->children[cIdx]->retrieveNode(r, depth);
}

/** Node retriever that ALWAYS returns the requested node, possibly without coefs.
  *
  * Recursive routine to find and return the node with a given NodeIndex. This
  * routine always returns the appropriate node, and will generate nodes that
  * does not exist. Recursion starts at at this node and ASSUMES the requested
  * node is in fact decending from this node. */
template<int D>
MWNode<D> *MWNode<D>::retrieveNode(const NodeIndex<D> &idx) {
    if (getScale() == idx.getScale()) { // we're done
        assert(getNodeIndex() == idx);
        return this;
    }
    assert(isAncestor(idx));
    lockNode();
    if (isLeafNode()) {
        genChildren();
        giveChildrenCoefs();
    }
    unlockNode();
    int cIdx = getChildIndex(idx);
    assert(this->children[cIdx] != 0);
    return this->children[cIdx]->retrieveNode(idx);
}

/** Test if a given coordinate is within the boundaries of the node. */
template<int D>
bool MWNode<D>::hasCoord(const double *r) const {
    double sFac = pow(2.0, -getScale());
    const int *l = getTranslation();
    //    println(1, "[" << r[0] << "," << r[1] << "," << r[2] << "]");
    //    println(1, "[" << l[0] << "," << l[1] << "," << l[2] << "]");
    //    println(1, *this);
    for (int d = 0; d < D; d++) {
        if (r[d] < sFac*l[d] or r[d] > sFac*(l[d] + 1)) {
            //            println(1, "false");
            return false;
        }
    }
    //    println(1, "true");
    return true;
}

/** Testing if nodes are compatible wrt NodeIndex and Tree (order, rootScale,
  * relPrec, etc). */
template<int D>
bool MWNode<D>::isCompatible(const MWNode<D> &node) {
    NOT_IMPLEMENTED_ABORT;
    //    if (nodeIndex != node.nodeIndex) {
    //        println(0, "nodeIndex mismatch" << std::endl);
    //        return false;
    //    }
    //    if (not this->tree->checkCompatible(*node.tree)) {
    //        println(0, "tree type mismatch" << std::endl);
    //        return false;
    //    }
    //    return true;
}

/** Test if the node is decending from a given NodeIndex, that is, if they have
  * overlapping support. */
template<int D>
bool MWNode<D>::isAncestor(const NodeIndex<D> &idx) const {
    int relScale = idx.getScale() - getScale();
    if (relScale < 0) {
        return false;
    }
    const int *l = getTranslation();
    for (int d = 0; d < D; d++) {
        int reqTransl = idx.getTranslation()[d] >> relScale;
        if (l[d] != reqTransl) {
            return false;
        }
    }
    return true;
}

template<int D>
bool MWNode<D>::isDecendant(const NodeIndex<D> &idx) const {
    NOT_IMPLEMENTED_ABORT;
}


//template<int D>
//void MWNode<D>::broadcastCoefs(int src, mpi::communicator *comm) {
    //NOT_IMPLEMENTED_ABORT;
    //#ifdef HAVE_MPI

    //    if (comm != 0) {
    //        comm = &node_group;
    //    }
    //    assert(this->isAllocated());
    //    double *data = coefs->data();
    //    mpi::broadcast(*comm, data, getNCoefs(), src);
    //    this->setHasCoefs();
    //    this->setFullRedundancy();
    //#endif
    //}

    //template<int D>
    //void MWNode<D>::sendCoefs(int dest, int tag,
    //                             mpi::communicator *comm) {
    //#ifdef HAVE_MPI

    //    if (comm != 0) {
    //        comm = &node_group;
    //    }
    //    const double *data = coefs->data();
    //    comm->send(dest, tag, data, getNCoefs());
    //    this->setRedundancy(dest);
    //#endif
//}

//template<int D>
//void MWNode<D>::assignDecendantTags(int rank) {
    //NOT_IMPLEMENTED_ABORT;
    //    for (int n = 0; n < getNChildren(); n++) {
    //        MWNode<D> &child = getMWChild(n);
    //        child.setRankId(rank);
    //        child.assignDecendantTags(rank);
    //    }
//}
template class MWNode<1>;
template class MWNode<2>;
template class MWNode<3>;
