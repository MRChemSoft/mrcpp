#include "TreeBuilder.h"
#include "TreeCalculator.h"
#include "TreeAdaptor.h"
#include "trees/MWTree.h"
#include "trees/MWNode.h"
#include "utils/Timer.h"
#include "utils/Printer.h"

using namespace std;

namespace mrcpp {

template<int D>
void TreeBuilder<D>::build(MWTree<D> &tree,
                           TreeCalculator<D> &calculator,
                           TreeAdaptor<D> &adaptor,
                           int maxIter) const {
    Timer calc_t(false), split_t(false), norm_t(false);
    println(10, " == Building tree");

    MWNodeVector *newVec = 0;
    MWNodeVector *workVec = calculator.getInitialWorkVector(tree);

    double sNorm = 0.0;
    double wNorm = 0.0;

    int iter = 0;
    while (workVec->size() > 0) {
        printout(10, "  -- #" << setw(3) << iter << ": Calculated ");
        printout(10, setw(6) << workVec->size() << " nodes ");
        calc_t.resume();
        calculator.calcNodeVector(*workVec);
        calc_t.stop();

        norm_t.resume();
        if (iter == 0) {
            sNorm = calcScalingNorm(*workVec);
        }
        wNorm += calcWaveletNorm(*workVec);

        if (sNorm < 0.0 or wNorm < 0.0) {
            tree.squareNorm = -1.0;
        } else {
            // approximate norm for thresholding only
            // exact norm is recomputed after mwTransform
            tree.squareNorm = sNorm + wNorm;
        }
        println(10, setw(24) << tree.squareNorm);
        norm_t.stop();

        split_t.resume();
        newVec = new MWNodeVector;
        if (iter >= maxIter and maxIter >= 0) workVec->clear();
        adaptor.splitNodeVector(*newVec, *workVec);
        split_t.stop();

        delete workVec;
        workVec = newVec;
        iter++;
    }
    tree.resetEndNodeTable();
    delete workVec;

    Printer::printSeparator(10, ' ');
    Printer::printTime(10, "Time calc", calc_t);
    Printer::printTime(10, "Time norm", norm_t);
    Printer::printTime(10, "Time split", split_t);
}

template<int D>
void TreeBuilder<D>::clear(MWTree<D> &tree, TreeCalculator<D> &calculator) const {
    println(10, " == Clearing tree");

    Timer clean_t;
    MWNodeVector nodeVec;
    tree.makeNodeTable(nodeVec);
    calculator.calcNodeVector(nodeVec);//clear all coefficients
    clean_t.stop();

    tree.clearSquareNorm();

    println(10, "  -- #  1: Cleared      " << setw(6) << nodeVec.size() << " nodes");
    Printer::printSeparator(10, ' ');
    Printer::printTime(10, "Time clean", clean_t);
    Printer::printSeparator(10, ' ');
}

template<int D>
int TreeBuilder<D>::split(MWTree<D> &tree, TreeAdaptor<D> &adaptor, bool passCoefs) const {
    println(10, " == Refining tree");

    Timer split_t;
    MWNodeVector newVec;
    MWNodeVector *workVec = tree.copyEndNodeTable();
    adaptor.splitNodeVector(newVec, *workVec);
    if (passCoefs) {
        for (int i = 0; i < workVec->size(); i++) {
            MWNode<D> &node = *(*workVec)[i];
            if (node.isBranchNode()) {
                node.giveChildrenCoefs(true);
            }
        }
    }
    delete workVec;
    tree.resetEndNodeTable();
    split_t.stop();

    printout(10, "  -- #  0: Split        ");
    printout(10, setw(6) << newVec.size() << " nodes\n");

    Printer::printSeparator(10, ' ');
    Printer::printTime(10, "Time split", split_t);
    Printer::printSeparator(10, ' ');

    return newVec.size();
}

template<int D>
void TreeBuilder<D>::calc(MWTree<D> &tree, TreeCalculator<D> &calculator) const {
    println(10, " == Calculating tree");

    Timer calc_t;
    MWNodeVector *workVec = calculator.getInitialWorkVector(tree);
    calculator.calcNodeVector(*workVec);
    printout(10, "  -- #" << setw(3) << 0 << ": Calculated ");
    printout(10, setw(6) << workVec->size() << " nodes ");
    delete workVec;
    calc_t.stop();

    tree.calcSquareNorm();

    Printer::printSeparator(10, ' ');
    Printer::printTime(10, "Time calc", calc_t);
}

template<int D>
double TreeBuilder<D>::calcScalingNorm(const MWNodeVector &vec) const {
    double sNorm = 0.0;
    for (int i = 0; i < vec.size(); i++) {
        const MWNode<D> &node = *vec[i];
        sNorm += node.getScalingNorm();
    }
    return sNorm;
}

template<int D>
double TreeBuilder<D>::calcWaveletNorm(const MWNodeVector &vec) const {
    double wNorm = 0.0;
    for (int i = 0; i < vec.size(); i++) {
        const MWNode<D> &node = *vec[i];
        wNorm += node.getWaveletNorm();
    }
    return wNorm;
}

template class TreeBuilder<1>;
template class TreeBuilder<2>;
template class TreeBuilder<3>;

} // namespace mrcpp
