#include "TreeBuilder.h"
#include "TreeCalculator.h"
#include "TreeAdaptor.h"
#include "MWTree.h"
#include "MWNode.h"

using namespace std;

template<int D>
TreeBuilder<D>::TreeBuilder(const MultiResolutionAnalysis<D> &mra, int iter)
        : adaptor(0),
          calculator(0),
          MRA(mra),
          maxIter(iter) {
}

template<int D>
TreeBuilder<D>::~TreeBuilder() {
    if (this->adaptor != 0) MSG_ERROR("Adaptor not deallocated");
    if (this->calculator != 0) MSG_ERROR("Calculator not deallocated");
}

template<int D>
void TreeBuilder<D>::clearAdaptor() {
    if (this->adaptor != 0) {
        delete this->adaptor;
        this->adaptor = 0;
    }
}

template<int D>
void TreeBuilder<D>::clearCalculator() {
    if (this->calculator != 0) {
        delete this->calculator;
        this->calculator = 0;
    }
}

template<int D>
void TreeBuilder<D>::build(MWTree<D> &tree) {
    if (this->calculator == 0) MSG_ERROR("Calculator not initialized");
    println(10, " == Building tree");

    NodeIndexSet *splitSet = 0;
    MWNodeVector *splitVec = 0;
    MWNodeVector *workVec = tree.copyEndNodeTable();
    MWNodeVector *endVec = tree.getEndNodeTable();
    endVec->clear();

    int iter = 0;
    bool computesCoefs = this->calculator->computesCoefs();
    if (not computesCoefs) tree.clearSquareNorm();
    while (workVec->size() > 0) {
        printout(10, "  -- #" << setw(3) << iter << ": Calculated   ");
        workVec = clearForeignNodes(workVec);
        this->calculator->calcNodeVector(*workVec);//set all coefficients
        if (computesCoefs) tree.calcSquareNorm(workVec);
        if (maxIterReached(iter) or this->adaptor == 0) break;
        splitVec = this->adaptor->splitNodeVector(*workVec, endVec);
        splitSet = getNodeIndexSet(*splitVec);
        broadcast_index_list<D>(*splitSet);
        workVec->clear();
        tree.splitNodes(*splitSet, workVec);//allocate new nodes
        delete splitVec;
        delete splitSet;
        iter++;
    }
    delete workVec;
    tree.resetEndNodeTable();
    if (computesCoefs) tree.calcSquareNorm();
}

template<int D>
MWNodeVector* TreeBuilder<D>::clearForeignNodes(MWNodeVector *oldVec) const {
    MWNodeVector *newVec = new MWNodeVector;
    for (int i = 0; i < oldVec->size(); i++) {
        MWNode<D> *node = (*oldVec)[i];
        if (node == 0) {
            continue;
        }
        if (not node->isForeign()) {
            newVec->push_back(node);
        }
    }
    delete oldVec;
    return newVec;
}

template<int D>
NodeIndexSet* TreeBuilder<D>::getNodeIndexSet(const MWNodeVector &nodeVec) const {
    NodeIndexSet *idxSet = new NodeIndexSet;
    for (int i = 0; i < nodeVec.size(); i++) {
        const NodeIndex<D> &idx = nodeVec[i]->getNodeIndex();
        idxSet->insert(&idx);
    }
    return idxSet;
}

template class TreeBuilder<1>;
template class TreeBuilder<2>;
template class TreeBuilder<3>;
