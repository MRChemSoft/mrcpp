#include "GridCleaner.h"
#include "MWTree.h"
#include "Timer.h"

using namespace std;

template<int D>
int GridCleaner<D>::clean(MWTree<D> &tree) const {
    Timer clean_t, split_t;
    if (this->calculator == 0) MSG_ERROR("Calculator not initialized");
    if (this->adaptor == 0) MSG_ERROR("Adaptor not initialized");
    println(10, " == Clearing tree");

    split_t.restart();
    MWNodeVector *newVec = new MWNodeVector;
    MWNodeVector *workVec = tree.copyEndNodeTable();
    this->adaptor->splitNodeVector(*newVec, *workVec);
    int nSplit = newVec->size();
    delete workVec;
    delete newVec;
    split_t.stop();

    printout(10, "  -- #  0: Split        ");
    printout(10, setw(6) << nSplit << " nodes\n");

    clean_t.restart();
    MWNodeVector nodeVec;
    tree.makeNodeTable(nodeVec);
    int nClear = nodeVec.size();
    this->calculator->calcNodeVector(nodeVec);//clear all coefficients
    clean_t.stop();

    tree.resetEndNodeTable();
    tree.clearSquareNorm();

    println(10, "  -- #  1: Cleared      " << setw(6) << nClear << " nodes");
    println(10, "");
    println(10, "Time split          " << split_t);
    println(10, "Time clean          " << clean_t);
    println(10, endl);

    return nSplit;
}

template class GridCleaner<1>;
template class GridCleaner<2>;
template class GridCleaner<3>;
