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

#include "TreeCalculator.h"

namespace mrcpp {

class ABGVCalculator final : public TreeCalculator<2> {
public:
    ABGVCalculator(const ScalingBasis &basis, double a, double b);

private:
    const double A; ///< Left boundary conditions, ref. Alpert et al.
    const double B; ///< Right boundary conditions, ref. Alpert et al.
    Eigen::MatrixXd K;
    Eigen::VectorXd valueZero;
    Eigen::VectorXd valueOne;

    void calcNode(MWNode<2> &node) override;

    void calcKMatrix(const ScalingBasis &basis);
    void calcValueVectors(const ScalingBasis &basis);
};

} // namespace mrcpp
