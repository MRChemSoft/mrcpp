#include <Eigen/Core>

#include "SerialTree.h"
#include "MWTree.h"
#include "TelePrompter.h"
#include "MathUtils.h"

using namespace Eigen;
using namespace std;

/** Overwrite all pointers defined in the tree.
  * Necessary after sending the tree 
  * could be optimized. Should reset other counters? (GenNodes...) */
template<int D>
void SerialTree<D>::rewritePointers(int nChunks){
    NOT_IMPLEMENTED_ABORT;
    /*
  int depthMax = 100;
  MWNode<D>* stack[depthMax*8];
  int slen = 0, counter = 0;

  this->nGenNodes = 0;
  this->nGenNodesCoeff = -1;
  this->nLooseNodesCoeff = 0;

  //reinitialize stacks
  for (int i = 0; i < this->maxNodes; i++) {
        this->nodeStackStatus[i] = 0;
  }

  for (int i = 0; i < this->maxGenNodes; i++) {
        this->genNodeStackStatus[i] = 0;//0=unoccupied
  }
  this->genNodeStackStatus[this->maxGenNodes] = -1;//=unavailable

  for (int i = 0; i < this->maxLooseNodesCoeff; i++) {
    this->looseCoeffStackStatus[i] = 0;//0=unoccupied
  }
  this->looseCoeffStackStatus[this->maxLooseNodesCoeff]=-1;//-1=unavailable

  this->getTree()->nNodes = 0;
  this->getTree()->nodesAtDepth.clear();
  this->getTree()->squareNorm = 0.0;

  for(int ichunk = 0 ; ichunk < nChunks; ichunk++){
    for(int inode = 0 ; inode < this->maxNodesPerChunk; inode++){
      ProjectedNode<D>* node = (this->nodeChunks[ichunk]) + inode;
      if (node->SNodeIx >= 0) {
	  //node is part of tree, should be processed
	  this->getTree()->incrementNodeCount(node->getScale());
	  if (node->isEndNode()) this->getTree()->squareNorm += node->getSquareNorm();
	  
	  //normally (intel) the virtual table does not change, but we overwrite anyway
	  *(char**)(Node) = this->cvptr_ProjectedNode;
	  
	  Node->tree = this->getTree();

	  //"adress" of coefs is the same as node, but in another array
	  Node->coefs = this->NodeCoeffChunks[ichunk]+ inode*this->sizeNodeCoeff;
	  
	  //adress of parent and children must be corrected
	  //can be on a different chunks
	  if(Node->parentSNodeIx>=0){
	    int n_ichunk = Node->parentSNodeIx/this->maxNodesPerChunk;
	    int n_inode = Node->parentSNodeIx%this->maxNodesPerChunk;
	    Node->parent = this->NodeChunks[n_ichunk] + n_inode;
	  }else{Node->parent = 0;}
	    
	  for (int i = 0; i < Node->getNChildren(); i++) {
	    int n_ichunk = (Node->childSNodeIx+i)/this->maxNodesPerChunk;
	    int n_inode = (Node->childSNodeIx+i)%this->maxNodesPerChunk;
	    Node->children[i] = this->NodeChunks[n_ichunk] + n_inode;
	  }
	  this->NodeStackStatus[Node->SNodeIx] = 1;//occupied
#ifdef OPENMP
	  omp_init_lock(&(Node->node_lock));
#endif
	}

    }
  }

  //update other MWTree data
  FunctionTree<D>* Tree = static_cast<FunctionTree<D>*> (this->mwTree_p);

  NodeBox<D> &rBox = Tree->getRootBox();
  MWNode<D> **roots = rBox.getNodes();

  for (int rIdx = 0; rIdx < rBox.size(); rIdx++) {
    roots[rIdx] = (this->NodeChunks[0]) + rIdx;//adress of roots are at start of NodeChunks[0] array
  }

  this->getTree()->resetEndNodeTable();

*/
}

/** Make children scaling coefficients from parent
 * Other node info are not used/set
 * coeff_in are not modified.
 * The output is written directly into the 8 children scaling coefficients. 
 * NB: ASSUMES that the children coefficients are separated by Children_Stride!
 */
template<int D>
void SerialTree<D>::S_mwTransform(double* coeff_in, double* coeff_out, bool readOnlyScaling, int stride, bool b_overwrite) {
    int operation = Reconstruction;
    int kp1 = this->getTree()->getKp1();
    int kp1_d = this->getTree()->getKp1_d();
    int tDim = (1<<D);
    int kp1_dm1 = MathUtils::ipow(kp1, D - 1);
    const MWFilter &filter = this->getTree()->getMRA().getFilter();
    double overwrite = 0.0;
    double *tmp;
    double tmpcoeff[kp1_d*tDim];
    double tmpcoeff2[kp1_d*tDim];
    int ftlim=tDim;
    int ftlim2=tDim;
    int ftlim3=tDim;
    if(readOnlyScaling){
        ftlim=1;
        ftlim2=2;
        ftlim3=4;
        //NB: Careful: tmpcoeff tmpcoeff2 are not initialized to zero
        //must not read these unitialized values!
    }

    overwrite = 0.0;
    int i = 0;
    int mask = 1;
    for (int gt = 0; gt < tDim; gt++) {
        double *out = tmpcoeff + gt * kp1_d;
        for (int ft = 0; ft < ftlim; ft++) {
            // Operate in direction i only if the bits along other
            // directions are identical. The bit of the direction we
            // operate on determines the appropriate filter/operator
            if ((gt | mask) == (ft | mask)) {
	        double *in = coeff_in + ft * kp1_d;
	        int filter_index = 2 * ((gt >> i) & 1) + ((ft >> i) & 1);
	        const MatrixXd &oper = filter.getSubFilter(filter_index, operation);

	        MathUtils::applyFilter(out, in, oper, kp1, kp1_dm1, overwrite);
	        overwrite = 1.0;
            }
        }
        overwrite = 0.0;
    }
    if (D>1) {
        i++;
        mask = 2;//1 << i;
        for (int gt = 0; gt < tDim; gt++) {
            double *out = tmpcoeff2 + gt * kp1_d;
            for (int ft = 0; ft < ftlim2; ft++) {
                // Operate in direction i only if the bits along other
                // directions are identical. The bit of the direction we
                // operate on determines the appropriate filter/operator
	        if ((gt | mask) == (ft | mask)) {
	            double *in = tmpcoeff + ft * kp1_d;
	            int filter_index = 2 * ((gt >> i) & 1) + ((ft >> i) & 1);
	            const MatrixXd &oper = filter.getSubFilter(filter_index, operation);
	  
	            MathUtils::applyFilter(out, in, oper, kp1, kp1_dm1, overwrite);
	            overwrite = 1.0;
	        }
            }
            overwrite = 0.0;
        }
    }
    if (D>2) {
        overwrite = 1.0;
        if(b_overwrite) overwrite = 0.0;
        i++;
        mask = 4;//1 << i;
        for (int gt = 0; gt < tDim; gt++) {
            double *out = coeff_out + gt * stride;//write right into children
            for (int ft = 0; ft < ftlim3; ft++) {
	        // Operate in direction i only if the bits along other
	        // directions are identical. The bit of the direction we
	        // operate on determines the appropriate filter/operator
                if ((gt | mask) == (ft | mask)) {
	            double *in = tmpcoeff2 + ft * kp1_d;
	            int filter_index = 2 * ((gt >> i) & 1) + ((ft >> i) & 1);
	            const MatrixXd &oper = filter.getSubFilter(filter_index, operation);

	            MathUtils::applyFilter(out, in, oper, kp1, kp1_dm1, overwrite);
	            overwrite = 1.0;
                }
            }
            overwrite = 1.0;
            if(b_overwrite) overwrite = 0.0;
        }
    }

    if (D>3) MSG_FATAL("D>3 NOT IMPLEMENTED for S_mwtransform");

    if (D<3) {
        double *out;
        if(D==1)out=tmpcoeff;
        if(D==2)out=tmpcoeff2;
        if(b_overwrite){
            for (int j = 0; j < tDim; j++){ 
	        for (int i = 0; i < kp1_d; i++){ 
	            coeff_out[i+j*stride] = out[i+j*kp1_d];
	        }
            }
        }else{
            for (int j = 0; j < tDim; j++){ 
	        for (int i = 0; i < kp1_d; i++){ 
	            coeff_out[i+j*stride]+=out[i+j*kp1_d];
	        }
            }
        }
    }
}

/** Make parent from children scaling coefficients
 * Other node info are not used/set
 * coeff_in are not modified.
 * The output is read directly from the 8 children scaling coefficients. 
 * NB: ASSUMES that the children coefficients are separated by Children_Stride!
 */
template<int D>
void SerialTree<D>::S_mwTransformBack(double* coeff_in, double* coeff_out, int stride) {
    NOT_IMPLEMENTED_ABORT;
}

template<>
void SerialTree<3>::S_mwTransformBack(double* coeff_in, double* coeff_out, int stride) {
  int operation = Compression;
  int kp1 = this->getTree()->getKp1();
  int kp1_d = this->getTree()->getKp1_d();
  int tDim = 8;
  int kp1_dm1 = MathUtils::ipow(kp1, 2);
  const MWFilter &filter = this->getTree()->getMRA().getFilter();
  double overwrite = 0.0;
  double tmpcoeff[kp1_d*tDim];

  int ftlim = tDim;
  int ftlim2 = tDim;
  int ftlim3 = tDim;

  int i = 0;
  int mask = 1;
  for (int gt = 0; gt < tDim; gt++) {
        double *out = coeff_out + gt * kp1_d;
        for (int ft = 0; ft < ftlim; ft++) {
            // Operate in direction i only if the bits along other
            // directions are identical. The bit of the direction we
            // operate on determines the appropriate filter/operator
            if ((gt | mask) == (ft | mask)) {
	        double *in = coeff_in + ft * stride;
	        int filter_index = 2 * ((gt >> i) & 1) + ((ft >> i) & 1);
	        const MatrixXd &oper = filter.getSubFilter(filter_index, operation);

	        MathUtils::applyFilter(out, in, oper, kp1, kp1_dm1, overwrite);
	        overwrite = 1.0;
            }
        }
        overwrite = 0.0;
    }
    i++;
    mask = 2;//1 << i;
    for (int gt = 0; gt < tDim; gt++) {
        double *out = tmpcoeff + gt * kp1_d;
        for (int ft = 0; ft < ftlim2; ft++) {
            // Operate in direction i only if the bits along other
            // directions are identical. The bit of the direction we
            // operate on determines the appropriate filter/operator
            if ((gt | mask) == (ft | mask)) {
	        double *in = coeff_out + ft * kp1_d;
	        int filter_index = 2 * ((gt >> i) & 1) + ((ft >> i) & 1);
	        const MatrixXd &oper = filter.getSubFilter(filter_index, operation);

	        MathUtils::applyFilter(out, in, oper, kp1, kp1_dm1, overwrite);
	        overwrite = 1.0;
            }
        }
        overwrite = 0.0;
    }
    i++;
    mask = 4;//1 << i;
    for (int gt = 0; gt < tDim; gt++) {
        double *out = coeff_out + gt * kp1_d;
        //double *out = coeff_out + gt * N_coeff;
        for (int ft = 0; ft < ftlim3; ft++) {
            // Operate in direction i only if the bits along other
            // directions are identical. The bit of the direction we
            // operate on determines the appropriate filter/operator
            if ((gt | mask) == (ft | mask)) {
	        double *in = tmpcoeff + ft * kp1_d;
	        int filter_index = 2 * ((gt >> i) & 1) + ((ft >> i) & 1);
	        const MatrixXd &oper = filter.getSubFilter(filter_index, operation);

	        MathUtils::applyFilter(out, in, oper, kp1, kp1_dm1, overwrite);
	        overwrite = 1.0;
            }
        }
        overwrite = 0.0;
    }
}

template class SerialTree<1>;
template class SerialTree<2>;
template class SerialTree<3>;
