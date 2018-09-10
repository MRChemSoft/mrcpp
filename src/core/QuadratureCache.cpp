/*
 *
 *
 *  \date Jul 26, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */

#include "QuadratureCache.h"
#include "utils/Printer.h"

namespace mrcpp {

QuadratureCache::QuadratureCache() {
    this->A = 0.0;
    this->B = 1.0;
    this->intervals = 1;
}

QuadratureCache::~QuadratureCache() {
}

void QuadratureCache::load(int k) {
    SET_CACHE_LOCK();
    if (not hasId(k)) {
        GaussQuadrature *gp = new GaussQuadrature(k, this->A, this->B, this->intervals);
        int memo = 2 * k * sizeof(double);
        ObjectCache<GaussQuadrature>::load(k, gp, memo);
    }
    UNSET_CACHE_LOCK();
}

GaussQuadrature &QuadratureCache::get(int k) {
    if (not hasId(k)) {
        load(k);
    }
    return ObjectCache<GaussQuadrature>::get(k);
}

void QuadratureCache::setBounds(double a, double b) {
    if (fabs(this->A - a) < MachineZero and fabs(this->B - b) < MachineZero) {
        return;
    }
    if (a >= b) {
        MSG_ERROR("Invalid Gauss interval, a > b.");
    }
    this->A = a;
    this->B = b;
    for (int i = 0; i < getNObjs(); i++) {
        if (hasId(i)) {
            ObjectCache<GaussQuadrature>::get(i).setBounds(a, b);
        }
    }
}

void QuadratureCache::setIntervals(int ivals) {
    if (ivals == this->intervals) {
        return;
    }
    if (this->intervals < 1) {
        MSG_ERROR("Invalid number of intervals, intervals < 1");
    }
    for (int i = 0; i < getNObjs(); i++) {
        if (hasId(i)) {
            ObjectCache<GaussQuadrature>::get(i).setIntervals(ivals);
        }
    }
}

} // namespace mrcpp
