#include "PHOperator.h"
#include "treebuilders/TreeBuilder.h"
#include "treebuilders/PHCalculator.h"
#include "treebuilders/BandWidthAdaptor.h"
#include "utils/Printer.h"
#include "utils/Timer.h"

namespace mrcpp {

template<int D>
PHOperator<D>::PHOperator(const MultiResolutionAnalysis<D> &mra, int order)
        : DerivativeOperator<D>(mra) {
    initializeOperator(order);
}

template<int D>
void PHOperator<D>::initializeOperator(int order) {
    int bw = 1; // Operator bandwidth
    int max_scale = this->oper_mra.getMaxScale();
    const ScalingBasis &basis = this->oper_mra.getScalingBasis();

    TreeBuilder<2> builder;
    PHCalculator calculator(basis, order);
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

template class PHOperator<1>;
template class PHOperator<2>;
template class PHOperator<3>;

} // namespace mrcpp
