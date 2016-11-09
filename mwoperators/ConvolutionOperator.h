#ifndef CONVOLUTIONOPERATOR_H
#define CONVOLUTIONOPERATOR_H

#include "MWOperator.h"
#include "FunctionTreeVector.h"

template<int D>
class ConvolutionOperator : public MWOperator {
public:
    ConvolutionOperator(const MultiResolutionAnalysis<D> &mra, double pr);
    virtual ~ConvolutionOperator();

protected:
    double prec;
    FunctionTreeVector<1> kernel_exp;
    MultiResolutionAnalysis<1> kern_mra;

    void initializeOperator(GreensKernel &greens_kernel);
    void clearKernel();

    double calcMinDistance(const MultiResolutionAnalysis<D> &MRA, double epsilon) const;
    double calcMaxDistance(const MultiResolutionAnalysis<D> &MRA) const;
};

#endif // CONVOLUTIONOPERATOR_H
