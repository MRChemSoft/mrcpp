#ifndef DERIVATIVECALCULATOR_H
#define DERIVATIVECALCULATOR_H

#include "TreeCalculator.h"

class DerivativeCalculator : public TreeCalculator<2> {
public:
    DerivativeCalculator(const ScalingBasis &basis, double a, double b);
    virtual ~DerivativeCalculator() { }

protected:
    const double A;	 ///< Left boundary conditions, ref. Alpert et al.
    const double B;	 ///< Right boundary conditions, ref. Alpert et al.
    Eigen::MatrixXd K;
    Eigen::VectorXd valueZero;
    Eigen::VectorXd valueOne;

    virtual void calcNode(MWNode<2> &node);

    void calcKMatrix(const ScalingBasis &basis);
    void calcValueVectors(const ScalingBasis &basis);
};

#endif // DERIVATIVECALCULATOR_H
