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
#include "LegendrePoly.h"

class LegendreBasis : public ScalingBasis {
public:
    LegendreBasis(int k)
            : ScalingBasis(k, Legendre) {
        initScalingBasis();
    }
    virtual ~LegendreBasis() { }

    void initScalingBasis() {
        for (int k = 0; k < getScalingOrder() + 1; k++) {
            LegendrePoly L_k(k, 2.0, 1.0);
            L_k *= sqrt(2.0 * k + 1.0); // exact normalization
            this->funcs.push_back(L_k);
        }
    }
};

#endif // LEGENDREBASIS_H
