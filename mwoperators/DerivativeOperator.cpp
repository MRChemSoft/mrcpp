#include "DerivativeOperator.h"
#include "DerivativeGenerator.h"

template<int D>
DerivativeOperator<D>::DerivativeOperator(int dir,
                                          const MultiResolutionAnalysis<D> &mra,
                                          double a,
                                          double b)
        : MWOperator<D>(mra, -1.0),
          A(a),
          B(b) {
    if (this->A > MachineZero) NEEDS_TESTING;
    if (this->B > MachineZero) NEEDS_TESTING;
    initializeOperator();
    this->setApplyDir(dir);
}

template<int D>
DerivativeOperator<D>::~DerivativeOperator() {
    this->clearOperator();
}

template<int D>
void DerivativeOperator<D>::initializeOperator() {
    MultiResolutionAnalysis<2> *oper_mra = this->getOperatorMRA();
    DerivativeGenerator DG(*oper_mra);

    OperatorTree *oper_comp = DG(this->A, this->B);
    this->oper.push_back(oper_comp);

    delete oper_mra;
}

template class DerivativeOperator<1>;
template class DerivativeOperator<2>;
template class DerivativeOperator<3>;
