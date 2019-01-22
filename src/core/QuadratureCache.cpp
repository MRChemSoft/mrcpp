/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2019 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
 *
 * This file is part of MRCPP.
 *
 * MRCPP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRCPP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRCPP.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRCPP, see:
 * <https://mrcpp.readthedocs.io/>
 */

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
        auto *gp = new GaussQuadrature(k, this->A, this->B, this->intervals);
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
    if (std::abs(this->A - a) < MachineZero and std::abs(this->B - b) < MachineZero) {
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
