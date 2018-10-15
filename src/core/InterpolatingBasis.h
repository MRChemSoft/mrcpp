#pragma once

#include "ScalingBasis.h"

namespace mrcpp {

class InterpolatingBasis final : public ScalingBasis {
public:
    InterpolatingBasis(int k) : ScalingBasis(k, Interpol) {
        initScalingBasis();
        calcQuadratureValues();
        calcCVMaps();
    }

    void initScalingBasis();
    void calcQuadratureValues();
    void calcCVMaps();
};

}
