#include "ConvolutionOperator.h"
#include "CrossCorrelationCalculator.h"
#include "OperatorAdaptor.h"
#include "GridGenerator.h"
#include "MWProjector.h"
#include "OperatorTree.h"
#include "GreensKernel.h"
#include "Gaussian.h"
#include "BandWidth.h"
#include "MathUtils.h"

using namespace std;
using namespace Eigen;

template<int D>
ConvolutionOperator<D>::ConvolutionOperator(const MultiResolutionAnalysis<D> &mra, double pr)
    : prec(pr),
      oper_mra(mra.getOperatorMRA()),
      kern_mra(mra.getKernelMRA()) {
}

template<int D>
ConvolutionOperator<D>::~ConvolutionOperator() {
    this->clearOperator();
    this->clearKernel();
}

template<int D>
void ConvolutionOperator<D>::initializeOperator(GreensKernel &greens_kernel) {
    int max_scale = this->oper_mra.getMaxScale();
    GridGenerator<1> G(max_scale);
    MWProjector<1> Q(this->prec/10.0, max_scale);

    TreeBuilder<2> builder;
    OperatorAdaptor adaptor(this->prec, max_scale);

    for (int i = 0; i < greens_kernel.size(); i++) {
        Gaussian<1> &k_func = *greens_kernel[i];
        FunctionTree<1> *k_tree = new FunctionTree<1>(this->kern_mra, MaxAllocNodes1D);
        G(*k_tree, k_func); //Generate empty grid to hold narrow Gaussian
        Q(*k_tree, k_func); //Project Gaussian starting from the empty grid
        CrossCorrelationCalculator calculator(*k_tree);

        OperatorTree *o_tree = new OperatorTree(this->oper_mra, this->prec);
        builder.build(*o_tree, calculator, adaptor, -1); //Expand 1D kernel into 2D operator

        Timer trans_t;
        o_tree->mwTransform(BottomUp);
        o_tree->calcSquareNorm();
        o_tree->setupOperNodeCache();
        trans_t.stop();

        println(10, "Time transform      " << trans_t);
        println(10, std::endl);

        this->kern_exp.push_back(k_tree);
        this->oper_exp.push_back(o_tree);
    }
}

template<int D>
void ConvolutionOperator<D>::clearOperator() {
    for (int i = 0; i < this->oper_exp.size(); i++) {
        if (this->oper_exp[i] != 0) delete this->oper_exp[i];
    }
    this->oper_exp.clear();
}

template<int D>
void ConvolutionOperator<D>::clearKernel() {
    for (int i = 0; i < this->kern_exp.size(); i++) {
        if (this->kern_exp[i] != 0) delete this->kern_exp[i];
    }
    this->kern_exp.clear();
}

template<int D>
double ConvolutionOperator<D>::calcMinDistance(const MultiResolutionAnalysis<D> &MRA, double epsilon) const {
    int maxScale = MRA.getMaxScale();
    return sqrt(epsilon * pow(2.0, -maxScale));
}

template<int D>
double ConvolutionOperator<D>::calcMaxDistance(const MultiResolutionAnalysis<D> &MRA) const {
    const double *lb = MRA.getWorldBox().getLowerBounds();
    const double *ub = MRA.getWorldBox().getUpperBounds();
    return MathUtils::calcDistance(D, lb, ub);
}

template<int D>
OperatorTree& ConvolutionOperator<D>::getComponent(int i) {
    if (this->oper_exp[i] == 0) MSG_ERROR("Invalid component");
    if (i < 0 or i >= this->oper_exp.size()) MSG_ERROR("Out of bounds");
    return *this->oper_exp[i];
}

template<int D>
const OperatorTree& ConvolutionOperator<D>::getComponent(int i) const {
    if (this->oper_exp[i] == 0) MSG_ERROR("Invalid component");
    if (i < 0 or i >= this->oper_exp.size()) MSG_ERROR("Out of bounds");
    return *this->oper_exp[i];
}

template<int D>
int ConvolutionOperator<D>::getMaxBandWidth(int depth) const {
    int maxWidth = -1;
    if (depth < 0 ) {
        maxWidth = this->bandMax.maxCoeff();
    } else if (depth < this->bandMax.size() ) {
        maxWidth = this->bandMax(depth);
    }
    return maxWidth;
}

template<int D>
void ConvolutionOperator<D>::calcBandWidths(double prec) {
    int maxDepth = 0;
    // First compute BandWidths and find depth of the deepest component
    for (unsigned int i = 0; i < this->oper_exp.size(); i++) {
        OperatorTree &oTree = *this->oper_exp[i];
        oTree.calcBandWidth(prec);
        const BandWidth &bw = oTree.getBandWidth();
        int depth = bw.getDepth();
        if (depth > maxDepth) {
            maxDepth = depth;
        }
    }
    this->bandMax = VectorXi(maxDepth + 1);
    this->bandMax.setConstant(-1);
    // Find the largest effective bandwidth at each scale
    for (unsigned int i = 0; i < this->oper_exp.size(); i++) {
        const OperatorTree &oTree = *this->oper_exp[i];
        const BandWidth &bw = oTree.getBandWidth();
        for (int n = 0; n <= bw.getDepth(); n++) { // scale loop
            for (int j = 0; j < 4; j++) { //component loop
                int w = bw.getWidth(n, j);
                if (w > this->bandMax(n)) {
                    this->bandMax(n) = w;
                }
            }
        }
    }
    println(20, "  Maximum bandwidths:\n" << this->bandMax << std::endl);
}

template<int D>
void ConvolutionOperator<D>::clearBandWidths() {
    for (unsigned int i = 0; i < this->oper_exp.size(); i++) {
        this->oper_exp[i]->clearBandWidth();
    }
}

template class ConvolutionOperator<1>;
template class ConvolutionOperator<2>;
template class ConvolutionOperator<3>;
