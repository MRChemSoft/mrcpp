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

#include <Eigen/Core>

#include "GaussQuadrature.h"
#include "ObjectCache.h"

namespace mrcpp {

#define getQuadratureCache(X) QuadratureCache &X = QuadratureCache::getInstance()

class QuadratureCache final : public ObjectCache<GaussQuadrature> {
public:
    static QuadratureCache &getInstance() {
        static QuadratureCache theQuadratureCache;
        return theQuadratureCache;
    }

    void load(int order);
    GaussQuadrature &get(int order);

    const Eigen::VectorXd &getRoots(int i) { return get(i).getRoots(); }
    const Eigen::VectorXd &getWeights(int i) { return get(i).getWeights(); }

    void setIntervals(int i);
    void setBounds(double a, double b);

    int getIntervals() const { return this->intervals; }
    double getUpperBound() const { return this->B; }
    double getLowerBound() const { return this->A; }

private:
    double A;
    double B;
    int intervals;

    QuadratureCache();
    ~QuadratureCache();

    QuadratureCache(QuadratureCache const &qc) = delete;
    QuadratureCache &operator=(QuadratureCache const &qc) = delete;
};

} // namespace mrcpp
