#ifndef GRIDCLEANER_H
#define GRIDCLEANER_H

#include "TreeBuilder.h"

template<int D>
class GridCleaner : public TreeBuilder<D> {
public:
    GridCleaner(const MultiResolutionAnalysis<D> &mra)
            : TreeBuilder<D>(mra, -1) {
        this->calculator = new TreeCalculator<D>();
    }
    
    GridCleaner(const MultiResolutionAnalysis<D> &mra,
                const TreeAdaptor<D> &a)
            : TreeBuilder<D>(mra, -1) {
        this->calculator = new TreeCalculator<D>();
        this->adaptor = a.copy();
    }

    virtual ~GridCleaner() {
        this->clearCalculator();
        this->clearAdaptor();
    }

    int operator()(FunctionTree<D> &out) {
        return clean(out);
    }

protected:

    int clean(MWTree<D> &tree) {
        if (this->calculator == 0) MSG_ERROR("Calculator not initialized");
        println(10, " == Clearing tree");

        int nSplit = 0;
        if (this->adaptor != 0) {
            MWNodeVector *workVec = tree.copyEndNodeTable();
            workVec = this->clearForeignNodes(workVec);
            MWNodeVector *splitVec = this->adaptor->splitNodeVector(*workVec);
            NodeIndexSet *splitSet = this->getNodeIndexSet(*splitVec);
            nSplit = splitVec->size();

            broadcast_index_list<D>(*splitSet);
            tree.splitNodes(*splitSet);//allocate new nodes

            delete workVec;
            delete splitVec;
            delete splitSet;
        }

        printout(10, "  -- #  0: Split        ");
        printout(10, std::setw(6) << nSplit << " nodes\n");
        printout(10, "  -- #  1: Cleared      ");

        MWNodeVector nodeVec;
        tree.makeNodeTable(nodeVec);
        this->calculator->calcNodeVector(nodeVec);//clear all coefficients

        tree.resetEndNodeTable();
        tree.clearSquareNorm();
        return nSplit;
    }
};

#endif // GRIDCLEANER_H
