#ifndef CONVOLUTIONOPERATOR_H
#define CONVOLUTIONOPERATOR_H

#include "MWOperator.h"
#include "FunctionTreeVector.h"

template<int D>
class ConvolutionOperator : public MWOperator<D> {
public:
    ConvolutionOperator(const MultiResolutionAnalysis<D> &mra,
                        double apply, double build, int dir);
    virtual ~ConvolutionOperator();

protected:
    double build_prec;
    FunctionTreeVector<1> kernel;

    void initializeOperator(GreensKernel &greens_kernel);
    void clearKernel();

    double calcMinDistance(double epsilon) const;
    double calcMaxDistance() const;

    MultiResolutionAnalysis<1> *getKernelMRA() const;
};

#endif // CONVOLUTIONOPERATOR_H
