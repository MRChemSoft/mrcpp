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

#include "ScalingBasis.h"

class LegendreBasis : public ScalingBasis {
public:
    LegendreBasis(int k)
            : ScalingBasis(k, Legendre) {
        initScalingBasis();
        calcQuadratureValues();
    }
    virtual ~LegendreBasis() { }

    void initScalingBasis();
    void calcQuadratureValues();
};

#endif // LEGENDREBASIS_H
