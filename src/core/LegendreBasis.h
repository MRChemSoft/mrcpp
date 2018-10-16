#pragma once

#pragma GCC system_header
#include <Eigen/Dense>

#include "ScalingBasis.h"

namespace mrcpp {

class LegendreBasis final : public ScalingBasis {
public:
    LegendreBasis(int k) : ScalingBasis(k, Legendre) {
        initScalingBasis();
        calcQuadratureValues();
        calcCVMaps();
    }

    void initScalingBasis();
    void calcQuadratureValues();
    void calcCVMaps();
};

}
