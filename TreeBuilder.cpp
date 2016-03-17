#include "TreeBuilder.h"
#include "TreeCalculator.h"
#include "TreeAdaptor.h"
#include "MWTree.h"
#include "MWNode.h"
#include "Timer.h"

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
void TreeBuilder<D>::build(MWTree<D> &tree) const {
    Timer calc_t, split_t;
    if (this->calculator == 0) MSG_ERROR("Calculator not initialized");
    if (this->adaptor == 0) MSG_ERROR("Adaptor not initialized");
    println(10, " == Building tree");

    MWNodeVector *newVec = 0;
    MWNodeVector *workVec = tree.copyEndNodeTable();

    double endNorm = 0.0;
    double workNorm = 0.0;

    int iter = 0;
    while (workVec->size() > 0) {
        printout(10, "  -- #" << setw(3) << iter++ << ": Calculated   ");
        printout(10, setw(6) << workVec->size() << " nodes\n");

        calc_t.restart();
        workNorm = this->calculator->calcNodeVector(*workVec);
        calc_t.stop();

        tree.squareNorm = max(endNorm + workNorm, -1.0);

        split_t.restart();
        newVec = new MWNodeVector;
        endNorm += this->adaptor->splitNodeVector(*newVec, *workVec);
        split_t.stop();

        delete workVec;
        workVec = newVec;
        if (maxIterReached(iter)) workVec->clear();
    }
    tree.resetEndNodeTable();
    delete workVec;

    println(10, "");
    println(10, "Time calc           " << calc_t);
    println(10, "Time split          " << split_t);
}

template<int D>
double TreeBuilder<D>::calcTreeNorm(MWNodeVector &workVec) const {
    double workNorm = 0.0;
    for (int i = 0; i < workVec.size(); i++) {
        MWNode<D> &node = *workVec[i];
        if (node.isEndNode()) {
            workNorm += node.getSquareNorm();
        }
    }
    return workNorm;
}

template class TreeBuilder<1>;
template class TreeBuilder<2>;
template class TreeBuilder<3>;
