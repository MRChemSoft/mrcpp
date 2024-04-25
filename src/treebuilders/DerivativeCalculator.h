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
#include "operators/OperatorStatistics.h"

namespace mrcpp {

template <int D, typename T> class DerivativeCalculator final : public TreeCalculator<D, T> {
public:
  DerivativeCalculator(int dir, DerivativeOperator<D> &o, FunctionTree<D, T> &f);
    ~DerivativeCalculator() override;

    MWNodeVector<D, T> *getInitialWorkVector(MWTree<D, T> &tree) const override;
    void calcNode(MWNode<D, T> &fNode, MWNode<D, T> &gNode);

private:
    int applyDir;
    FunctionTree<D, T> *fTree;
    DerivativeOperator<D> *oper;

    std::vector<Timer> band_t;
    std::vector<Timer> calc_t;
    std::vector<Timer> norm_t;
    OperatorStatistics<D, T> operStat;

    MWNodeVector<D, T> makeOperBand(const MWNode<D, T> &gNode, std::vector<NodeIndex<D>> &idx_band);

    void initTimers();
    void clearTimers();
    void printTimers() const;

    void calcNode(MWNode<D, T> &node) override;
    void postProcess() override {
        printTimers();
        clearTimers();
        initTimers();
    }

    void applyOperator(OperatorState<D, T> &os);
    void applyOperator_bw0(OperatorState<D, T> &os);
    void tensorApplyOperComp(OperatorState<D, T> &os);
};

} // namespace mrcpp
