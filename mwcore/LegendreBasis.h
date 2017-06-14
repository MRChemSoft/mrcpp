/*
 *
 *
 *  \date June 2, 2010
 *  \author Stig Rune Jensen \n
 *          CTCC, University of Tromsø
 *
 */

#ifndef LEGENDREBASIS_H
#define LEGENDREBASIS_H

#pragma GCC system_header
#include <Eigen/Dense>

#include "ScalingBasis.h"

class LegendreBasis : public ScalingBasis {
public:
    LegendreBasis(int k) : ScalingBasis(k, Legendre) {
        initScalingBasis();
        calcQuadratureValues();
        calcCVMaps();
    }
    virtual ~LegendreBasis() { }

    void initScalingBasis();
    void calcQuadratureValues();
    void calcCVMaps();
};

#endif // LEGENDREBASIS_H
