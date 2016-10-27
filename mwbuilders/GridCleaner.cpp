#include "GridCleaner.h"
#include "MWTree.h"
#include "Timer.h"

using namespace std;

template<int D>
int GridCleaner<D>::clean(MWTree<D> &tree,
                          TreeCalculator<D> &calculator,
                          TreeAdaptor<D> &adaptor) const {
    println(10, " == Clearing tree");

    Timer split_t;
    MWNodeVector newVec;
    MWNodeVector *workVec = tree.copyEndNodeTable();
    adaptor.splitNodeVector(newVec, *workVec);
    int nSplit = newVec.size();
    delete workVec;
    split_t.stop();

    printout(10, "  -- #  0: Split        ");
    printout(10, setw(6) << nSplit << " nodes\n");

    Timer clean_t;
    MWNodeVector nodeVec;
    tree.makeNodeTable(nodeVec);
    int nClear = nodeVec.size();
    calculator.calcNodeVector(nodeVec);//clear all coefficients
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
