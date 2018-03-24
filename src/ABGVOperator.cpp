#include "ABGVOperator.h"
#include "ABGVCalculator.h"
#include "TreeBuilder.h"
#include "BandWidthAdaptor.h"
#include "OperatorTree.h"
#include "Printer.h"
#include "Timer.h"

using namespace mrcpp;

template<int D>
ABGVOperator<D>::ABGVOperator(const MultiResolutionAnalysis<D> &mra, double a, double b)
        : DerivativeOperator<D>(mra) {
    initializeOperator(a, b);
}

template<int D>
void ABGVOperator<D>::initializeOperator(double a, double b) {
    int bw = 0; // Operator bandwidth
    if (fabs(a) > MachineZero) bw = 1;
    if (fabs(b) > MachineZero) bw = 1;
    int max_scale = this->oper_mra.getMaxScale();
    const ScalingBasis &basis = this->oper_mra.getScalingBasis();

    TreeBuilder<2> builder;
    ABGVCalculator calculator(basis, a, b);
    BandWidthAdaptor adaptor(bw, max_scale);

    OperatorTree *o_tree = new OperatorTree(this->oper_mra, MachineZero);
    builder.build(*o_tree, calculator, adaptor, -1);

    Timer trans_t;
    o_tree->calcSquareNorm();
    o_tree->setupOperNodeCache();
    trans_t.stop();

    println(10, "Time transform      " << trans_t);
    println(10, std::endl);

    this->oper_exp.push_back(o_tree);
}

template class mrcpp::ABGVOperator<1>;
template class mrcpp::ABGVOperator<2>;
template class mrcpp::ABGVOperator<3>;
