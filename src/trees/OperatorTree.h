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

#pragma once

#include "MWTree.h"

namespace mrcpp {

class OperatorTree final : public MWTree<2> {
public:
    OperatorTree(const MultiResolutionAnalysis<2> &mra, double np);
    OperatorTree(const OperatorTree &tree) = delete;
    OperatorTree &operator=(const OperatorTree &tree) = delete;
    virtual ~OperatorTree();

    double getNormPrecision() const { return this->normPrec; }

    void calcBandWidth(double prec = -1.0);
    void clearBandWidth();

    void setupOperNodeCache();
    void clearOperNodeCache();

    BandWidth &getBandWidth() { return *this->bandWidth; }
    const BandWidth &getBandWidth() const { return *this->bandWidth; }

    OperatorNode &getNode(int n, int l) { return *nodePtrAccess[n][l]; }
    const OperatorNode &getNode(int n, int l) const { return *nodePtrAccess[n][l]; }

    void mwTransformDown(bool overwrite);
    void mwTransformUp();

protected:
    const double normPrec;
    BandWidth *bandWidth;
    OperatorNode ***nodePtrStore;  ///< Avoids tree lookups
    OperatorNode ***nodePtrAccess; ///< Center (l=0) of node list

    void getMaxTranslations(Eigen::VectorXi &maxTransl);

    std::ostream &print(std::ostream &o);
};

} // namespace mrcpp
