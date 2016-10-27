#include "DerivativeOperator.h"
#include "DerivativeGenerator.h"
#include "GridGenerator.h"
#include "MWAdder.h"

template<int D>
DerivativeOperator<D>::DerivativeOperator(const MultiResolutionAnalysis<D> &mra,
                                          double a, double b)
        : MWOperator<D>(mra),
          A(a),
          B(b) {
    if (this->A > MachineZero) NEEDS_TESTING;
    if (this->B > MachineZero) NEEDS_TESTING;
    initializeOperator();
}

template<int D>
DerivativeOperator<D>::~DerivativeOperator() {
    this->clearOperator();
}

template<int D>
void DerivativeOperator<D>::initializeOperator() {
    MultiResolutionAnalysis<2> *oper_mra = this->MRA.getOperatorMRA();
    DerivativeGenerator DG(this->MRA.getScalingBasis());

    OperatorTree *oper_comp = new OperatorTree(*oper_mra, MachineZero, MaxAllocOperNodes);
    DG(*oper_comp, this->A, this->B);
    this->oper_exp.push_back(oper_comp);

    delete oper_mra;
}

template<int D>
void DerivativeOperator<D>::grad(FunctionTreeVector<D> &out,
                                 FunctionTree<D> &inp) {
    NOT_IMPLEMENTED_ABORT;
    /*
    if (out.size() != 0) MSG_ERROR("Invalid input");

    GridGenerator<D> G;
    for (int d = 0; d < D; d++) {
        this->setApplyDir(d);
        FunctionTree<D> *out_d = new FunctionTree<D>(this->MRA);
        G(*out_d, inp);
        (*this)(*out_d, inp, 0);
        out.push_back(out_d);
    }
    */
}

template<int D>
void DerivativeOperator<D>::div(FunctionTree<D> &out,
                                FunctionTreeVector<D> &inp) {
    NOT_IMPLEMENTED_ABORT;
    /*
    if (inp.size() != D) MSG_ERROR("Invalid input");

    GridGenerator<D> G;
    MWAdder<D> add;

    FunctionTreeVector<D> vec;
    for (int d = 0; d < D; d++) {
        this->setApplyDir(d);
        FunctionTree<D> *tmp_d = new FunctionTree<D>(this->MRA);
        G(*tmp_d, *inp[d]);
        (*this)(*tmp_d, *inp[d], 0);
        vec.push_back(tmp_d);
    }

    G(out, vec);
    add(out, vec, 0);
    vec.clear(true);
    */
}

template class DerivativeOperator<1>;
template class DerivativeOperator<2>;
template class DerivativeOperator<3>;
