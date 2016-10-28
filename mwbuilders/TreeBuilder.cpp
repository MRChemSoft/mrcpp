#include "TreeBuilder.h"
#include "TreeCalculator.h"
#include "TreeAdaptor.h"
#include "MWTree.h"
#include "MWNode.h"
#include "Timer.h"

using namespace std;

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
        println(10, setw(25) << tree.squareNorm);
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

    println(10, "");
    println(10, "Time calc           " << calc_t);
    println(10, "Time norm           " << norm_t);
    println(10, "Time split          " << split_t);
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
