#ifndef DERIVATIVECONVOLUTION_H
#define DERIVATIVECONVOLUTION_H

#include "ConvolutionOperator.h"
#include "GridGenerator.h"
#include "DerivativeKernel.h"

template<int D>
class DerivativeConvolution : public ConvolutionOperator<D> {
public:
    DerivativeConvolution(int d,
                          const MultiResolutionAnalysis<D> &mra,
                          double apply = -1.0,
                          double build = -1.0)
            : ConvolutionOperator<D>(mra, apply, build) {
        double epsilon = this->build_prec/10.0;
        DerivativeKernel derivative_kernel(epsilon);
        this->initializeOperator(derivative_kernel);
        this->setApplyDir(d);
    }
    virtual ~DerivativeConvolution() { }

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
    FunctionTreeVector<3> curl(FunctionTreeVector<3> &inp) {
        NOT_IMPLEMENTED_ABORT;
    }
};

#endif // DERIVATIVECONVOLUTION_H
