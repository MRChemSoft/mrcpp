/**
 *  \date April 20, 2010
 *  CTCC, University of Troms
 *
 */

#include "MWTree.h"
#include "MultiResolutionAnalysis.h"
#include "HilbertIterator.h"
#include "MathUtils.h"

using namespace std;
using namespace Eigen;

/** MWTree constructor.
  * Creates an empty tree object. Node construction and assignment of most of
  * the parameters are done in derived classes. */
template<int D>
MWTree<D>::MWTree(const MultiResolutionAnalysis<D> &mra)
        : nThreads(omp_get_max_threads()),
          MRA(mra),
          rootBox(mra.getWorldBox()),
          order(mra.getOrder()),
          kp1_d(MathUtils::ipow(mra.getOrder() + 1, D)),
          squareNorm(-1.0),
          name("nn"),
          nNodes(0),
          serialTree_p(0) {
    this->nodesAtDepth.push_back(0);
    allocNodeCounters();
    allocWorkMemory();
    println(10, "new MWTree ");
#ifdef OPENMP
    omp_init_lock(&tree_lock);
#endif
}

/** MWTree constructor with SerialTree storage for nodes.
  * Creates an empty tree object. Node construction and assignment of most of
  * the parameters are done in derived classes. */
template<int D>
MWTree<D>::MWTree(const MultiResolutionAnalysis<D> &mra, int max_nodes)
        : nThreads(omp_get_max_threads()),
          MRA(mra),
          rootBox(mra.getWorldBox()),
          order(mra.getOrder()),
          kp1_d(MathUtils::ipow(mra.getOrder() + 1, D)),
          squareNorm(-1.0),
          name("nn"),
          nNodes(0) {
    this->nodesAtDepth.push_back(0);
    allocNodeCounters();
    allocWorkMemory();
    new SerialTree<D>(this, max_nodes);
    println(10, "new Serial MWTree ");
#ifdef OPENMP
    omp_init_lock(&tree_lock);
#endif
}

/** MWTree copy constructor.
  * Takes the parameters of the input tree, not it's data */
template<int D>
MWTree<D>::MWTree(const MWTree<D> &tree)
        : nThreads(omp_get_max_threads()),
          MRA(tree.MRA),
          rootBox(tree.rootBox),
          order(tree.order),
          kp1_d(tree.kp1_d),
          squareNorm(-1.0),
          name("nn"),
          nNodes(0) {
    this->nodesAtDepth.push_back(0);
    allocNodeCounters();

#ifdef OPENMP
    omp_init_lock(&tree_lock);
#endif
}

/** MWTree destructor. */
template<int D>
MWTree<D>::~MWTree() {
    println(10, "~MWTree");
    if(this->serialTree_p) {
      //SerialTree removes nodes
      println(10, "delete serialTree_p");
      delete this->serialTree_p;
    } else {
      //Has to remove nodes here?
      //MWNode<D> **roots = this->getRootBox().getNodes();
      //for (int i = 0; i < this->getRootBox().size(); i++) {
      //  //ProjectedNode<D> *node = static_cast<ProjectedNode<D> *>(roots[i]);
      //  roots[i]->~MWNode();
      //  roots[i] = 0;
      //}
    }
    this->endNodeTable.clear();
    if (this->nNodes != 0) {
      MSG_ERROR("Node count != 0 -> " << this->nNodes);
    }
    if (this->nodesAtDepth.size() != 1) {
      MSG_ERROR("Nodes at depth != 1 -> " << this->nodesAtDepth.size());
    }
    if (this->nodesAtDepth[0] != 0) {
      MSG_ERROR("Nodes at depth 0 != 0 -> " << this->nodesAtDepth[0]);
    }
    deleteNodeCounters();
<<<<<<< HEAD
    freeWorkMemory();
    
=======

>>>>>>> 543587dec7bf4b87c66e6cbb6d15a8ef42555c44
#ifdef OPENMP
    omp_destroy_lock(&tree_lock);
#endif

}

template<int D>
double MWTree<D>::estimateError(bool absPrec) {
    NOT_IMPLEMENTED_ABORT;
//    double error = 0.0;
//    for (int i = 0; i < this->getNEndNodes(); i++) {
//        MWNode<D> &node = getEndMWNode(i);
//        error += node.estimateError(absPrec);
//    }
//#ifdef HAVE_MPI
//    error = mpi::all_reduce(node_group, error, std::plus<double>());
//#endif
//    return error;
}

/** Calculate the squared norm of a function represented as a tree.
  *
  * Norm is calculated using endNodes only, but if your endNodeTable is
  * incomplete (e.g. within growTree), the missing nodes must be given in the
  * input work vector. Involves an MPI reduction operation. */
template<int D>
void MWTree<D>::calcSquareNorm() {
    double treeNorm = 0.0;
    for (int n = 0; n < this->getNEndNodes(); n++) {
        const MWNode<D> &node = getEndMWNode(n);
        assert(node.hasCoefs());
        treeNorm += node.getSquareNorm();
    }
    this->squareNorm = treeNorm;
}

/** Reduce the accuracy of the tree by deleting nodes
  * which have a higher precision than the requested precison.
  * By default, the relative precision of the tree is used. */
template<int D>
void MWTree<D>::crop(double thrs, bool absPrec) {
    NOT_IMPLEMENTED_ABORT;
}

//template<int D>
//void MWTree<D>::cropTree(double prec, bool absPrec) {
//    set<const NodeIndex<D> *, NodeIndexComp<D> > cropNodes;
//    for (int i = 0; i < this->rootBox.getNBoxes(); i++) {
//        MWNode<D> &rootNode = getRootMWNode(i);
//        if (this->isScattered()) {
//            rootNode.cropChildren(prec, &cropNodes);
//        } else {
//            rootNode.cropChildren(prec);
//        }
//    }
//    if (this->isScattered()) {
//        broadcastIndexList(cropNodes);
//        typename set<const NodeIndex<D> *,
//                NodeIndexComp<D> >::reverse_iterator it;
//        for (it = cropNodes.rbegin(); it != cropNodes.rend(); it++) {
//            MWNode<D> *node = this->findNode(**it);
//            if (node != 0) {
//                node->deleteChildren();
//            }
//        }
//    }
//    resetEndNodeTable();
//    this->squareNorm = calcSquareNorm();
//}

template<int D>
void MWTree<D>::mwTransform(int type, bool overwrite) {
    switch(type) {
    case TopDown:
        mwTransformDown(overwrite);
        break;
    case BottomUp:
        mwTransformUp(overwrite);
        break;
    default:
        MSG_FATAL("Invalid wavelet transform");
    }
}

/** Regenerate all s/d-coeffs by backtransformation, starting at the bottom and
  * thus purifying all coefficients. Option to overwrite or add up existing
  * coefficients of BranchNodes (can be used after operator application). */
template<int D>
void MWTree<D>::mwTransformUp(bool overwrite) {
    vector<MWNodeVector > nodeTable;
    this->makeNodeTable(nodeTable);
//#pragma omp parallel firstprivate(overwrite) shared(nodeTable)
//{
    int start = nodeTable.size() - 2;
    for (int n = start; n >= 0; n--) {
        int nNodes = nodeTable[n].size();
//#pragma omp for schedule(guided)
            for (int i = 0; i < nNodes; i++) {
                MWNode<D> &node = *nodeTable[n][i];
                if (node.isBranchNode()) {
                    node.reCompress(overwrite);
                }
            }
    }
//}
}

/** Regenerate all scaling coeffs by MW transformation of existing s/w-coeffs
  * on coarser scales, starting at the rootNodes. Option to overwrite or add up
  * existing scaling coefficients (can be used after operator application). */
template<int D>
void MWTree<D>::mwTransformDown(bool overwrite) {
    vector<MWNodeVector > nodeTable;
    makeNodeTable(nodeTable);
//#pragma omp parallel shared(nodeTable)
//        {
    for (int n = 0; n < nodeTable.size(); n++) {
        int n_nodes = nodeTable[n].size();
//#pragma omp for schedule(guided)
        for (int i = 0; i < n_nodes; i++) {
            MWNode<D> &node = *nodeTable[n][i];
            if (node.isBranchNode()) {
                node.giveChildrenCoefs(overwrite);
            }
        }
    }
//        }
}

/** Traverse tree and set all nodes to zero.
  *
  * Keeps the node structure of the tree, even though the zero function is
  * representable at depth zero. Use cropTree to remove unnecessary nodes.*/
template<int D>
void MWTree<D>::setZero() {
    HilbertIterator<D> it(this);
    while(it.next()) {
        MWNode<D> &node = it.getNode();
        node.zeroCoefs();
    }
    this->squareNorm = 0.0;
}

template<int D>
void MWTree<D>::allocNodeCounters() {
    this->nGenNodes = new int[this->nThreads];
    this->nAllocGenNodes = new int[this->nThreads];
    for (int i = 0; i < this->nThreads; i++) {
        this->nGenNodes[i] = 0;
        this->nAllocGenNodes[i] = 0;
    }
}

template<int D>
void MWTree<D>::deleteNodeCounters() {
    delete[] this->nGenNodes;
    delete[] this->nAllocGenNodes;
}

/** Increment node counters for non-GenNodes. This routine is not thread
  * safe, and must NEVER be called outside a critical region in parallel.
  * It's way. way too expensive to lock the tree, so don't even think
  * about it. */
template<int D>
void MWTree<D>::incrementNodeCount(int scale) {
    int depth = scale - getRootScale();
    assert(depth >= 0);
    int n = this->nodesAtDepth.size() - 1;
    if (depth > n) {
        for (int i = 0; i < depth - n; i++) {
            this->nodesAtDepth.push_back(0);
        }
    }
    this->nodesAtDepth[depth]++;
    this->nNodes++;
}

/** Decrement node counters for non-GenNodes. This routine is not thread
  * safe, and must NEVER be called outside a critical region in parallel.
  * It's way. way too expensive to lock the tree, so don't even think
  * about it. */
template<int D>
void MWTree<D>::decrementNodeCount(int scale) {
    int depth = scale - getRootScale();
    assert(depth >= 0);
    assert(depth < this->nodesAtDepth.size());
    this->nodesAtDepth[depth]--;
    assert(this->nodesAtDepth[depth] >= 0);
    if (this->nodesAtDepth[depth] == 0 and this->nodesAtDepth.size() > 1) {
        this->nodesAtDepth.pop_back();
    }
    this->nNodes--;
    assert(this->nNodes >= 0);
}

/** Update GenNode counts in a safe way. Since GenNodes are created on the
  * fly, we cannot control when to update the node counters without locking
  * the whole tree. Therefore GenNodes update thread-private counters, which
  * get merged with the correct global counters in xxxNodes[0]. This method
  * should be called outside of the parallel region for performance reasons. */
template<int D>
void MWTree<D>::updateGenNodeCounts() {
    lockTree();
    for (int i = 1; i < this->nThreads; i++) {
        this->nGenNodes[0] += this->nGenNodes[i];
        this->nAllocGenNodes[0] += this->nAllocGenNodes[i];
        this->nGenNodes[i] = 0;
        this->nAllocGenNodes[i] = 0;
    }
    assert(this->nGenNodes[0] >= 0);
    assert(this->nAllocGenNodes[0] >= 0);
    unlockTree();
}

/** Adds a GenNode to the count. */
template<int D>
void MWTree<D>::incrementGenNodeCount() {
    int n = omp_get_thread_num();
    assert(n >= 0);
    assert(n < this->nThreads);
    this->nGenNodes[n]++;
}

/** Removes a GenNode from the count. */
template<int D>
void MWTree<D>::decrementGenNodeCount() {
    int n = omp_get_thread_num();
    assert(n >= 0);
    assert(n < this->nThreads);
    this->nGenNodes[n]--;
}

/** Adds an allocated GenNode to the count. */
template<int D>
void MWTree<D>::incrementAllocGenNodeCount() {
    int n = omp_get_thread_num();
    assert(n >= 0);
    assert(n < this->nThreads);
    this->nAllocGenNodes[n]++;
}

/** Removes an allocated GenNode from the count. */
template<int D>
void MWTree<D>::decrementAllocGenNodeCount() {
    int n = omp_get_thread_num();
    assert(n >= 0);
    assert(n < this->nThreads);
    this->nAllocGenNodes[n]--;
}

/** Get Node count. */
template<int D>
int MWTree<D>::getNNodes(int depth) const {
    if (depth < 0) {
        return this->nNodes;
    }
    if (depth >= this->nodesAtDepth.size()) {
        return 0;
    }
    return this->nodesAtDepth[depth];
}

/** Get allocated GenNode count. Includes an OMP reduction operation. */
template<int D>
int MWTree<D>::getNAllocGenNodes() {
    updateGenNodeCounts();
    return this->nAllocGenNodes[0];
}

/** Get GenNode count. Includes an OMP reduction operation. */
template<int D>
int MWTree<D>::getNGenNodes() {
    updateGenNodeCounts();
    return this->nGenNodes[0];
}

/** Find and return the node with the given NodeIndex, const version.
  *
  * Recursive routine to find and return the node with a given NodeIndex.
  * This routine returns the appropriate ProjectedNode, or a NULL pointer if
  * the node does not exist, or if it is a GenNode. Recursion starts at the
  * appropriate rootNode. */
template<int D>
const MWNode<D>* MWTree<D>::findNode(const NodeIndex<D> &idx) const {
    const MWNode<D> &root = getRootBox().getNode(idx);
    assert(root.isAncestor(idx));
    return root.retrieveNodeNoGen(idx);
}

/** Find and return the node with the given NodeIndex.
  *
  * Recursive routine to find and return the node with a given NodeIndex.
  * This routine returns the appropriate ProjectedNode, or a NULL pointer if
  * the node does not exist, or if it is a GenNode. Recursion starts at the
  * appropriate rootNode. */
template<int D>
MWNode<D>* MWTree<D>::findNode(const NodeIndex<D> &idx) {
    MWNode<D> &root = this->rootBox.getNode(idx);
    assert(root.isAncestor(idx));
    return root.retrieveNodeNoGen(idx);
}

/** Find and return the node with the given NodeIndex.
  *
  * This routine ALWAYS returns the node you ask for, and will generate nodes
  * that does not exist. Recursion starts at the appropriate rootNode and
  * decends from this.*/
template<int D>
MWNode<D>& MWTree<D>::getNode(const NodeIndex<D> &idx) {
    MWNode<D> &root = getRootBox().getNode(idx);
    assert(root.isAncestor(idx));
    return *root.retrieveNode(idx);
}

/** Find and return the node with the given NodeIndex.
  *
  * This routine returns the ProjectedNode you ask for, or the EndNode on
  * the path to the requested node, and will never create or return GenNodes.
  * Recursion starts at the appropriate rootNode and decends from this. */
template<int D>
MWNode<D>& MWTree<D>::getNodeOrEndNode(const NodeIndex<D> &idx) {
    MWNode<D> &root = getRootBox().getNode(idx);
    assert(root.isAncestor(idx));
    return *root.retrieveNodeOrEndNode(idx);
}

/** Find and return the node with the given NodeIndex.
  *
  * This routine returns the ProjectedNode you ask for, or the EndNode on
  * the path to the requested node, and will never create or return GenNodes.
  * Recursion starts at the appropriate rootNode and decends from this. */
template<int D>
const MWNode<D>& MWTree<D>::getNodeOrEndNode(const NodeIndex<D> &idx) const {
    const MWNode<D> &root = getRootBox().getNode(idx);
    assert(root.isAncestor(idx));
    return *root.retrieveNodeOrEndNode(idx);
}

/** Find and return the node at a given depth that contains a given coordinate.
  *
  * This routine ALWAYS returns the node you ask for, and will generate nodes
  * that does not exist. Recursion starts at the appropriate rootNode and
  * decends from this. */
template<int D>
MWNode<D>& MWTree<D>::getNode(const double *r, int depth) {
    MWNode<D> &root = getRootBox().getNode(r);
    if (depth >= 0) {
        return *root.retrieveNode(r, depth);
    } else {
        return *root.retrieveNodeOrEndNode(r, depth);
    }
}

/** Find and return the node at a given depth that contains a given coordinate.
  *
  * This routine returns the ProjectedNode you ask for, or the EndNode on
  * the path to the requested node, and will never create or return GenNodes.
  * Recursion starts at the appropriate rootNode and decends from this. */
template<int D>
MWNode<D>& MWTree<D>::getNodeOrEndNode(const double *r, int depth) {
    MWNode<D> &root = getRootBox().getNode(r);
    return *root.retrieveNodeOrEndNode(r, depth);
}

/** Find and return the node at a given depth that contains a given coordinate.
  *
  * This routine returns the ProjectedNode you ask for, or the EndNode on
  * the path to the requested node, and will never create or return GenNodes.
  * Recursion starts at the appropriate rootNode and decends from this. */
template<int D>
const MWNode<D>& MWTree<D>::getNodeOrEndNode(const double *r, int depth) const {
    const MWNode<D> &root = getRootBox().getNode(r);
    return *root.retrieveNodeOrEndNode(r, depth);
}

/** Traverse tree along the Hilbert path and find nodes of any rankId.
  * Returns one nodeVector for the whole tree. GenNodes disregarded. */
template<int D>
void MWTree<D>::makeNodeTable(MWNodeVector &nodeTable) {
    HilbertIterator<D> it(this);
    it.setReturnGenNodes(false);
    while (it.next()) {
        MWNode<D> &node = it.getNode();
        nodeTable.push_back(&node);
    }
}

/** Traverse tree along the Hilbert path and find nodes of any rankId.
  * Returns one nodeVector per scale. GenNodes disregarded. */
template<int D>
void MWTree<D>::makeNodeTable(std::vector<MWNodeVector > &nodeTable) {
    HilbertIterator<D> it(this);
    it.setReturnGenNodes(false);
    while (it.next()) {
        MWNode<D> &node = it.getNode();
        int depth = node.getDepth();
        if (depth + 1 > nodeTable.size()) { // Add one more element
            nodeTable.push_back(MWNodeVector());
        }
        nodeTable[depth].push_back(&node);
    }
}

template<int D>
MWNodeVector* MWTree<D>::copyEndNodeTable() {
    MWNodeVector *nVec = new MWNodeVector;
    for (int n = 0; n < getNEndNodes(); n++) {
        MWNode<D> &node = getEndMWNode(n);
        nVec->push_back(&node);
    }
    return nVec;
}

template<int D>
void MWTree<D>::resetEndNodeTable() {
    clearEndNodeTable();
    HilbertIterator<D> it(this);
    it.setReturnGenNodes(false);
    while (it.next()) {
        MWNode<D> &node = it.getNode();
        if (node.isEndNode()) {
            this->endNodeTable.push_back(&node);
        }
    }
}

template<int D>
int MWTree<D>::countBranchNodes(int depth) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
int MWTree<D>::countLeafNodes(int depth) {
    NOT_IMPLEMENTED_ABORT;
//    int nNodes = 0;
//    HilbertIterator<D> it(this);
//    while (it.next()) {
//        MWNode<D> &node = it.getNode();
//        if (node.getDepth() == depth or depth < 0) {
//            if (node.isLeafNode()) {
//                nNodes++;
//            }
//        }
//    }
//    return nNodes;
}

/** Traverse tree and count nodes. */
template<int D>
void MWTree<D>::RecountNodes() {
  int DepthMax = 100, slen = 0;
  MWNode<D>* stack[DepthMax*8];
  for (int rIdx = 0; rIdx < this->getRootBox().size(); rIdx++) {
   stack[slen++] = this->getRootBox().getNodes()[rIdx];
  }
  this->nNodes = 0;
  while (slen) {
    this->nNodes++;
    MWNode<D>* fpos = stack[--slen];
    for (int i = 0; i < fpos->getNChildren(); i++)stack[slen++] = fpos->children[i];
  }
}
/** Traverse tree and count nodes belonging to this rank. */
template<int D>
int MWTree<D>::countNodes(int depth) {
    NOT_IMPLEMENTED_ABORT;
//    HilbertIterator<D> it(this);
//    int count = 0;
//    while (it.next()) {
//        MWNode<D> &node = it.getNode();
//        if (node.isGenNode()) {
//            continue;
//        }
//        if (not node.isForeign()) {
//            count++;
//        }
//    }
//    return count;
}

/** Traverse tree and count nodes with allocated coefficients. */
template<int D>
int MWTree<D>::countAllocNodes(int depth) {
    NOT_IMPLEMENTED_ABORT;
//    HilbertIterator<D> it(this);
//    int count = 0;
//    while (it.next()) {
//        MWNode<D> &node = it.getNode();
//        if (node.isGenNode()) {
//            continue;
//        }
//        if (node.hasCoefs()) {
//            count++;
//        }
//    }
//    return count;
}

template<int D>
void MWTree<D>::deleteGenerated() {
    for (int n = 0; n < getNEndNodes(); n++) {
        getEndMWNode(n).deleteGenerated();
    }
}


/** Loop through endNodeTable and recursively clear all GenNode coefficients.
  * Includes a static cast of endNodes from MWNode to FunctionNode*/
template<int D>
void MWTree<D>::clearGenerated() {
    for (int i = 0; i < this->endNodeTable.size(); i++) {
        getEndMWNode(i).clearGenerated();
    }
}

template class MWTree<1>;
template class MWTree<2>;
template class MWTree<3>;
