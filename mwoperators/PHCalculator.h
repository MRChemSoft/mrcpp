#ifndef PHCALCULATOR_H
#define PHCALCULATOR_H

#include <Eigen/Core>

#include "TreeCalculator.h"

class PHCalculator : public TreeCalculator<2> {
public:
    PHCalculator(const ScalingBasis &basis);
    virtual ~PHCalculator() { }

protected:
    Eigen::MatrixXd S_m1;
    Eigen::MatrixXd S_0;
    Eigen::MatrixXd S_p1;

    virtual void calcNode(MWNode<2> &node);
    void readSMatrix(const ScalingBasis &basis);
};

#endif // PHCALCULATOR_H
