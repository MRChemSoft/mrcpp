#ifndef NUMDIFFCALCULATOR_H
#define NUMDIFFCALCULATOR_H

#include <Eigen/Core>

#include "TreeCalculator.h"

class NumDiffCalculator : public TreeCalculator<2> {
public:
    NumDiffCalculator(const ScalingBasis &basis);
    virtual ~NumDiffCalculator() { }

protected:
    Eigen::MatrixXd S_m1;
    Eigen::MatrixXd S_0;
    Eigen::MatrixXd S_p1;

    virtual void calcNode(MWNode<2> &node);
    void readSMatrix(const ScalingBasis &basis);
};

#endif // NUMDIFFCALCULATOR_H
