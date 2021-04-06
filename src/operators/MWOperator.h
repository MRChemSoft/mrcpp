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

#include <vector>

#include "trees/MultiResolutionAnalysis.h"
#include "trees/OperatorTree.h"

namespace mrcpp {

class MWOperator {
public:
    MWOperator(const MultiResolutionAnalysis<2> &mra)
            : oper_mra(mra) {}
    MWOperator(const MWOperator &oper) = delete;
    MWOperator &operator=(const MWOperator &oper) = delete;
    virtual ~MWOperator() { this->clear(true); }

    int size() const { return this->oper_exp.size(); }
    void push_back(OperatorTree *oper) { this->oper_exp.push_back(oper); }
    void clear(bool dealloc = false);

    int getMaxBandWidth(int depth = -1) const;
    const Eigen::VectorXi &getMaxBandWidths() const { return this->band_max; }

    void calcBandWidths(double prec);
    void clearBandWidths();

    OperatorTree &getComponent(int i);
    const OperatorTree &getComponent(int i) const;

    OperatorTree *operator[](int i) { return this->oper_exp[i]; }
    const OperatorTree *operator[](int i) const { return this->oper_exp[i]; }

protected:
    MultiResolutionAnalysis<2> oper_mra;
    OperatorTreeVector oper_exp;
    Eigen::VectorXi band_max;
};

} // namespace mrcpp
