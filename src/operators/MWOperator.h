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

template <int D>
class MWOperator {
public:
    MWOperator(const MultiResolutionAnalysis<D> &mra, int root, int reach);
    MWOperator(const MWOperator &oper) = delete;
    MWOperator &operator=(const MWOperator &oper) = delete;
    virtual ~MWOperator() = default;

    int size() const { return this->oper_exp.size(); }
    void push_back(std::unique_ptr<OperatorTree> oper) { this->oper_exp.push_back(std::move(oper)); }

    int getMaxBandWidth(int depth = -1) const;
    const std::vector<int> &getMaxBandWidths() const { return this->band_max; }

    void calcBandWidths(double prec);
    void clearBandWidths();

    int getOperatorRoot() const { return this->oper_root; }
    int getOperatorReach() const { return this->oper_reach; }

    OperatorTree &getComponent(int i);
    const OperatorTree &getComponent(int i) const;

    OperatorTree &operator[](int i) { return *this->oper_exp[i]; }
    const OperatorTree &operator[](int i) const { return *this->oper_exp[i]; }

protected:
    int oper_root;
    int oper_reach;
    MultiResolutionAnalysis<D> MRA;
    std::vector<std::unique_ptr<OperatorTree>> oper_exp;
    std::vector<int> band_max;

    MultiResolutionAnalysis<2> getOperatorMRA() const;
};

} // namespace mrcpp
