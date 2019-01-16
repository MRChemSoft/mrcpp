#pragma once

#pragma GCC system_header
#include <Eigen/Core>

#include "TreeCalculator.h"

namespace mrcpp {

class BSCalculator final : public TreeCalculator<2> {
public:
    explicit BSCalculator(const ScalingBasis &basis, int n);

private:
    const int diff_order;
    Eigen::MatrixXd S_m1;
    Eigen::MatrixXd S_0;
    Eigen::MatrixXd S_p1;

    void calcNode(MWNode<2> &node);
    void readSMatrix(const ScalingBasis &basis, char n);
};

}
