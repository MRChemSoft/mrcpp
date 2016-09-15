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
          tmpCoefs(0) {
    this->nodesAtDepth.push_back(0);
    allocNodeCounters();
    allocWorkMemory();
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
    allocWorkMemory();

#ifdef OPENMP
    omp_init_lock(&tree_lock);
#endif
}

/** MWTree destructor. */
template<int D>
MWTree<D>::~MWTree() {
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
    freeWorkMemory();

#ifdef OPENMP
    omp_destroy_lock(&tree_lock);
#endif
}

/** Allocate work memory of the tree, for use in mwTransform and the like. */
template<int D>
void MWTree<D>::allocWorkMemory() {
    int n_coefs = this->getTDim()*this->getKp1_d();
    this->tmpCoefs = new double*[this->nThreads];
    for (int i = 0; i < this->nThreads; i++) {
        this->tmpCoefs[i] = new double[n_coefs];
    }

    /*
    this->tmpCoefs = new Eigen::MatrixXd *[this->nThreads];
    this->tmpVector = new Eigen::VectorXd *[this->nThreads];
    this->tmpMWCoefs = new Eigen::VectorXd *[this->nThreads];
    for (int i = 0; i < this->nThreads; i++) {
        this->tmpCoefs[i] = new Eigen::MatrixXd(this->order + 1, D);
        this->tmpVector[i] = new Eigen::VectorXd(this->kp1_d);
        this->tmpMWCoefs[i] = new Eigen::VectorXd(this->kp1_d * (1 << D));
    }
    */
}

/** Deallocate work memory */
template<int D>
void MWTree<D>::freeWorkMemory() {
    if (this->tmpCoefs != 0) {
        for (int i = 0; i < this->nThreads; i++) {
            if (this->tmpCoefs[i] != 0) {
                delete[] this->tmpCoefs[i];
                this->tmpCoefs[i] = 0;
            }
        }
        delete[] this->tmpCoefs;
        this->tmpCoefs = 0;
    }
    /*
    for (int i = 0; i < this->nThreads; i++) {
        delete this->tmpCoefs[i];
        delete this->tmpVector[i];
        delete this->tmpMWCoefs[i];
    }
    delete[] this->tmpCoefs;
    delete[] this->tmpVector;
    delete[] this->tmpMWCoefs;
    */
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
#pragma omp parallel firstprivate(overwrite) shared(nodeTable)
{
    int start = nodeTable.size() - 2;
    for (int n = start; n >= 0; n--) {
        int nNodes = nodeTable[n].size();
#pragma omp for schedule(guided)
            for (int i = 0; i < nNodes; i++) {
                MWNode<D> &node = *nodeTable[n][i];
                if (node.isBranchNode()) {
                    node.reCompress(overwrite);
                }
            }
    }
}
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

/** Split nodes according to a list of NodeIndices.
  *
  * Given a list of NodeIndices to split, this routine creates the new children
  * nodes. The newly born (local) children nodes are collected in a MWNodeVector.
  * Children nodes are by default given the rank of their parent.*/
//template<int D>
//void MWTree<D>::splitNodes(const NodeIndexSet &idxSet, MWNodeVector *nVec) {
//    typename set<const NodeIndex<D> *>::iterator it;
//    for (it = idxSet.begin(); it != idxSet.end(); it++) {
//        MWNode<D> &node = getNode(**it);
//        node.createChildren();
//        if (nVec != 0) {
//            for (int i = 0; i < node.getNChildren(); i++) {
//                MWNode<D> *child = &node.getMWChild(i);
//                nVec->push_back(child);
//            }
//        }
//    }
//}

/** Increment node counters for non-GenNodes. This routine is not thread
  * safe, and must NEVER be called outside a critical region in parallel.
  * It's way. way too expensive to lock the tree, so don't even think
  * about it. */
template<int D>
void MWTree<D>::incrementNodeCount(int scale) {
    lockTree();
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
    unlockTree();
}

/** Decrement node counters for non-GenNodes. This routine is not thread
  * safe, and must NEVER be called outside a critical region in parallel.
  * It's way. way too expensive to lock the tree, so don't even think
  * about it. */
template<int D>
void MWTree<D>::decrementNodeCount(int scale) {
    lockTree();
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
    unlockTree();
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

/** Traverse nodeTable and find all nodes of different rankId. */
//template<int D>
//void MWTree<D>::findMissingNodes(MWNodeVector &nodeTable,
//                                 set<MWNode<D> *> &missing) {
//    NOT_IMPLEMENTED_ABORT;
//    for (unsigned int i = 0; i < nodeTable.size(); i++) {
//        MWNode<D> &node = *nodeTable[i];
//        if (not node.hasCoefs()) {
//            assert(node.isForeign());
//            missing.insert(&node);
//        }
//    }
//}

/** Traverse nodeTable and find all nodes with parent of different rankId. */
//template<int D>
//void MWTree<D>::findMissingParents(MWNodeVector &nodeTable,
//                                   set<MWNode<D> *> &missing) {
//    NOT_IMPLEMENTED_ABORT;
//    for (unsigned int i = 0; i < nodeTable.size(); i++) {
//        MWNode<D> &node = *nodeTable[i];
//        if (node.isRoot()) {
//            continue;
//        }
//        MWNode<D> &parent = node.getMWParent();
//        if (not parent.hasCoefs()) {
//            assert(this->getRankId() != parent.getRankId());
//            missing.insert(&parent);
//        }
//    }
//}

/** Traverse nodeTable and find all nodes with children of different rankId. */
//template<int D>
//void MWTree<D>::findMissingChildren(
//        MWNodeVector &nodeTable, set<MWNode<D> *> &missing) {
//    NOT_IMPLEMENTED_ABORT;
//    for (unsigned int i = 0; i < nodeTable.size(); i++) {
//        MWNode<D> &node = *nodeTable[i];
//        if (node.isEndNode()) {
//            continue;
//        }
//        for (int n = 0; n < node.getNChildren(); n++) {
//            MWNode<D> &child = node.getMWChild(n);
//            if (not child.hasCoefs()) {
//                assert(this->getRankId() != child.getRankId());
//                missing.insert(&child);
//            }
//        }
//    }
//}

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

/** Traverse tree along the Hilbert path and find nodes of local rankId.
  * Returns one nodeVector for the whole tree. GenNodes disregarded. */
//template<int D>
//void MWTree<D>::makeLocalNodeTable(MWNodeVector &nodeTable, bool common) {
//    NOT_IMPLEMENTED_ABORT;
//    HilbertIterator<D> it(this);
//    while (it.next()) {
//        MWNode<D> &node = it.getNode();
//        if (node.isGenNode()) {
//            continue;
//        }
//        if (node.isLocal() or (node.isCommon() and common)) {
//            nodeTable.push_back(&node);
//        }
//    }
//}

/** Traverse tree along the Hilbert path and find nodes of local rankId.
  * Returns one nodeVector per scale. GenNodes disregarded. */
//template<int D>
//void MWTree<D>::makeLocalNodeTable(std::vector<MWNodeVector > &nodeTable, bool common) {
//    HilbertIterator<D> it(this);
//    while (it.next()) {
//        MWNode<D> &node = it.getNode();
//        if (node.isGenNode()) {
//            continue;
//        }
//        int depth = node.getDepth();
//        if (depth + 1 > nodeTable.size()) { // Add one more element
//            nodeTable.push_back(MWNodeVector());
//        }
//        if (node.isLocal() or (node.isCommon() and common)) {
//            nodeTable[depth].push_back(&node);
//        }
//    }
//}

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

/** Print the number of nodes, sorted by depth and MPI rank. */
//template<int D>
//void MWTree<D>::printNodeRankCount() {
//    NOT_IMPLEMENTED_ABORT;
//    int nHosts = node_group.size();
//    int mDepth = getDepth();
//    int count[mDepth + 1][nHosts+1];
//    for (int i = 0; i < mDepth + 1; i++) {
//        for (int j = 0; j < nHosts+1; j++) {
//            count[i][j] = 0;
//        }
//    }
//    HilbertIterator<D> it(this);
//    while(it.next()) {
//        MWNode<D> &node = it.getNode();
//        if (node.isGenNode()) {
//            continue;
//        }
//        int depth = node.getDepth();
//        int rank = node.getRankId();
//        if (rank >= 0) {
//            count[depth][rank+1]++;
//            count[mDepth][rank+1]++;
//        } else {
//            count[depth][0]++;
//            count[mDepth][0]++;
//        }
//    }

//    println(0, endl);
//    printout(0, "   |");
//    for (int j = -1; j < nHosts; j++) {
//        printout(0, setw(6) << j);
//    }
//    println(0, "|" << endl);
//    for (int i = 0; i < mDepth + 1; i++) {
//        if (i == mDepth) {
//            printout(0, endl);
//        }
//        printout(0, setw(3) << i << "|");
//        for (int j = 0; j < nHosts+1; j++) {
//            printout(0, setw(6) << count[i][j]);
//        }
//        println(0, "|");
//    }
//    println(0, endl);
//}

//template<int D>
//void MWTree<D>::distributeNodes(int depth) {
//    NOT_IMPLEMENTED_ABORT;
//    MWNodeVector nodeTable;
//    HilbertIterator<D> it(this);
//    it.setReturnGenNodes(false);
//    it.setMaxDepth(depth);
//    while (it.next()) {
//        MWNode<D> &node = it.getNode();
//        if (node.isEndNode() or node.getDepth() == depth) {
//            nodeTable.push_back(&node);
//        }
//    }
//    distributeNodeTags(nodeTable);
//    tagDecendants(nodeTable);
//}

/** Tag each node with the rank who owns it. */
//template<int D>
//void MWTree<D>::distributeNodeTags(MWNodeVector &nodeList) {
//    NOT_IMPLEMENTED_ABORT;
//    int start, end;
//    int nNodes = nodeList.size();
//    int nHosts = node_group.size();
//    for (int k = 0; k < nHosts; k++) {
//        get_locale_index_range(k, nNodes, start, end);
//        for (int i = start; i < end; i++) {
//            MWNode<D> &node = *nodeList[i];
//            node.setRankId(k);
//        }
//    }
//}

/** Set the rank id (node tag) for all nodes in list. */
//template<int D>
//void MWTree<D>::tagNodes(MWNodeVector &nodeList, int rank) {
//    NOT_IMPLEMENTED_ABORT;
//    for (int i = 0; i < nodeList.size(); i++) {
//        MWNode<D> &node = *nodeList[i];
//        node.setRankId(rank);
//    }
//}

//template<int D>
//void MWTree<D>::tagDecendants(MWNodeVector &nodeList) {
//    NOT_IMPLEMENTED_ABORT;
//    int nNodes = nodeList.size();
//    for (int i = 0; i < nNodes; i++) {
//        MWNode<D> &node = *nodeList[i];
//        node.assignDecendantTags(node.getRankId());
//    }
//}

/** Traverse tree and remove nodes of foreign rank.
  * Option to keep all endNodes. */
//template<int D>
//void MWTree<D>::deleteForeign(bool keepEndNodes) {
//    NOT_IMPLEMENTED_ABORT;
//    if (not this->isScattered()) {
//        return;
//    }
//    HilbertIterator<D> it(this);

//    while (it.next()) {
//        MWNode<D> &node = it.getNode();
//        if (keepEndNodes and node.isEndNode()) {
//            continue;
//        }
//        if (node.isForeign()) {
//            node.clearCoefs();
//            node.clearNorms();
//        }
//        node.clearRedundancy();
//    }
//}

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

//template<int D>
//void MWTree<D>::checkGridOverlap(MWTree<D> &tree) {
//    NOT_IMPLEMENTED_ABORT;
//    int overlapA = 0;
//    int overlapB = 0;
//    int onlyA = 0;
//    int onlyB = 0;

//    HilbertIterator<D> itA(this);
//    itA.setReturnGenNodes(false);
//    while (itA.next()) {
//        const NodeIndex<D> &idx = itA.getNode().getNodeIndex();
//        MWNode<D> *nodeB = tree.findNode(idx);
//        if (nodeB != 0) {
//            overlapA++;
//        } else {
//            onlyA++;
//        }
//    }

//    HilbertIterator<D> itB(&tree);
//    itB.setReturnGenNodes(false);
//    while (itB.next()) {
//        const NodeIndex<D> &idx = itB.getNode().getNodeIndex();
//        MWNode<D> *nodeA = this->findNode(idx);
//        if (nodeA != 0) {
//            overlapB++;
//        } else {
//            onlyB++;
//        }
//    }

//    int nodesA = this->getNNodes();
//    int nodesB = tree.getNNodes();

//    if (overlapA != overlapB) {
//        MSG_WARN("Something went wrong, overlaps do not match.");
//    }
//    if (nodesA != (overlapA + onlyA)) {
//        MSG_WARN("Something went wrong, overlaps do not match.");
//    }
//    if (nodesB != (overlapB + onlyB)) {
//        MSG_WARN("Something went wrong, overlaps do not match.");
//    }

//    printout(0, "Overlapping nodes: ");
//    println(0, setw(8) << onlyA << setw(8) << overlapA << setw(8) << onlyB);
//}

//template<int D>
//void MWTree<D>::checkRankOverlap(MWTree<D> &tree) {
//    NOT_IMPLEMENTED_ABORT;
//    MatrixXi rankDiff = MatrixXi::Zero(2,11);

//    for (int i = 0; i < 11; i++) {
//        rankDiff(0,i) = -5 + i;
//    }

//    HilbertIterator<D> itA(this);
//    itA.setReturnGenNodes(false);
//    while (itA.next()) {
//        MWNode<D> &nodeA = itA.getNode();
//        MWNode<D> &nodeB = tree.getNodeNoGen(nodeA.getNodeIndex());
//        int rankA = nodeA.getRankId();
//        int rankB = nodeB.getRankId();
//        int diff = rankA - rankB;
//        if (diff <= -5) {
//            diff = 0;
//        } else if (diff >= 5) {
//            diff = 10;
//        } else {
//            diff += 5;
//        }
//        rankDiff(1, diff)++;
//    }

//    HilbertIterator<D> itB(&tree);
//    itB.setReturnGenNodes(false);
//    while (itB.next()) {
//        MWNode<D> &nodeB = itB.getNode();
//        MWNode<D> &nodeA = tree.getNodeNoGen(nodeB.getNodeIndex());
//        if (nodeA.getScale() == nodeB.getScale()) {
//            continue;
//        }
//        int rankA = nodeA.getRankId();
//        int rankB = nodeB.getRankId();
//        int diff = rankA - rankB;
//        if (diff <= -5) {
//            diff = 0;
//        } else if (diff >= 5) {
//            diff = 10;
//        } else {
//            diff += 5;
//        }
//        rankDiff(1, diff)++;
//    }
//    println(0, rankDiff.row(1));
//}

template class MWTree<1>;
template class MWTree<2>;
template class MWTree<3>;
