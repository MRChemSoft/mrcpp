/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2021 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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

#pragma once

#include "ObjectCache.h"

namespace mrcpp {

#define getLegendreScalingCache(X) ScalingCache<LegendreBasis> &X = ScalingCache<LegendreBasis>::getInstance()
#define getInterpolatingScalingCache(X)                                                                                \
    ScalingCache<InterpolatingBasis> &X = ScalingCache<InterpolatingBasis>::getInstance()

template <class P> class ScalingCache final : public ObjectCache<P> {
public:
    static ScalingCache &getInstance() {
        static ScalingCache theScalingCache;
        return theScalingCache;
    }
    void load(int order) {
        MRCPP_SET_OMP_LOCK();
        if (not this->hasId(order)) {
            P *f = new P(order);
            int memo = 2 * SQUARE(order + 1) * sizeof(double); // approx
            ObjectCache<P>::load(order, f, memo);
        }
        MRCPP_UNSET_OMP_LOCK();
    }

    P &get(int order) {
        if (not this->hasId(order)) { load(order); }
        return ObjectCache<P>::get(order);
    }

private:
    ScalingCache() {}
    ScalingCache(const ScalingCache<P> &sc) = delete;
    ScalingCache<P> &operator=(const ScalingCache<P> &sc) = delete;
};

} // namespace mrcpp
