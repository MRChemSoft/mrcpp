#ifndef DERIVATIVEOPERATOR_H
#define DERIVATIVEOPERATOR_H

#include "MWOperator.h"
#include "MWAdder.h"
#include "GridGenerator.h"

template<int D>
class DerivativeOperator : public MWOperator<D> {
public:
    DerivativeOperator(int dir,
                       const MultiResolutionAnalysis<D> &mra,
                       double a = 0.5,
                       double b = 0.5);
    virtual ~DerivativeOperator();

    FunctionTreeVector<D> grad(FunctionTree<D> &inp) {
        GridGenerator<D> G(this->MRA);
        FunctionTreeVector<D> out;
        for (int d = 0; d < D; d++) {
            this->setApplyDir(d);
            FunctionTree<D> *out_d = G(inp);
            (*this)(*out_d, inp, 0);
            out.push_back(out_d);
        }
        return out;
    }
    FunctionTree<D>* div(FunctionTreeVector<D> &inp) {
        if (inp.size() != D) MSG_ERROR("Invalid dimension");

        MWAdder<D> add(this->MRA);
        GridGenerator<D> G(this->MRA);

        FunctionTreeVector<D> vec;
        for (int d = 0; d < D; d++) {
            this->setApplyDir(d);
            FunctionTree<D> *tmp_d = G(inp);
            (*this)(*tmp_d, *inp[d], 0);
            vec.push_back(tmp_d);
        }

        FunctionTree<D> *out = G(inp);
        add(*out, vec, 0);
        vec.clear(true);

        return out;
    }
protected:
    const double A;
    const double B;

    void initializeOperator();
};

#endif // DERIVATIVEOPERATOR_H
