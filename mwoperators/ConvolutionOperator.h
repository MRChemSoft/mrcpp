#ifndef CONVOLUTIONOPERATOR_H
#define CONVOLUTIONOPERATOR_H

#include "MWOperator.h"
#include "FunctionTreeVector.h"

template<int D>
class ConvolutionOperator : public MWOperator<D> {
public:
    ConvolutionOperator(const MultiResolutionAnalysis<D> &mra, double pr)
        : MWOperator<D>(mra), prec(pr) { }
    virtual ~ConvolutionOperator();

protected:
    double prec;
    FunctionTreeVector<1> kernel_exp;

    void initializeOperator(GreensKernel &greens_kernel);
    void clearKernel();

    double calcMinDistance(double epsilon) const;
    double calcMaxDistance() const;
};

#endif // CONVOLUTIONOPERATOR_H
