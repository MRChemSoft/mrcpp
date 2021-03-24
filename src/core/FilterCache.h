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

/*
 *
 * \breif FilterCache is a static class taking care of loading and
 * unloading MultiWavelet filters, and their tensor counter parts.
 *
 * All data in FilterCache is static, and thus shared amongst all
 * instance objects. The type of filter, Legendre or Interpolating is
 * determined by a template variable so that both types of filters can
 * co-exist.
 *
 */

#pragma once

#include "MWFilter.h"
#include "ObjectCache.h"

#include <iostream>
#include <string>

namespace mrcpp {

#define getFilterCache(T, X) FilterCache<T> &X = FilterCache<T>::getInstance()
#define getLegendreFilterCache(X) FilterCache<Legendre> &X = FilterCache<Legendre>::getInstance()
#define getInterpolatingFilterCache(X) FilterCache<Interpol> &X = FilterCache<Interpol>::getInstance()

/** This class is an abstract base class for the various filter caches.
 * It's needed in order to be able to use the actual filter caches
 * without reference to the actual filter types */
class BaseFilterCache : public ObjectCache<MWFilter> {
public:
    void load(int order) override = 0;
    MWFilter &get(int order) override = 0;
    virtual const Eigen::MatrixXd &getFilterMatrix(int order) = 0;
};

template <int T> class FilterCache final : public BaseFilterCache {
public:
    static FilterCache &getInstance() {
        static FilterCache theFilterCache;
        return theFilterCache;
    }

    void load(int order) override;
    MWFilter &get(int order) override;
    const Eigen::MatrixXd &getFilterMatrix(int order) override;

protected:
    int type;

private:
    FilterCache();
    FilterCache(FilterCache<T> const &fc) = delete;
    FilterCache &operator=(FilterCache<T> const &fc) = delete;
};

} // namespace mrcpp
