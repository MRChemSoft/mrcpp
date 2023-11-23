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

template <int D> class DerivativeCalculator final : public TreeCalculator<D> {
public:
    DerivativeCalculator(int dir, DerivativeOperator<D> &o, FunctionTree<D> &f);
    ~DerivativeCalculator() override;

    MWNodeVector<D> *getInitialWorkVector(MWTree<D> &tree) const override;
    void calcNode(MWNode<D> &fNode, MWNode<D> &gNode);

private:
    int applyDir;
    FunctionTree<D> *fTree;
    DerivativeOperator<D> *oper;

    std::vector<Timer> band_t;
    std::vector<Timer> calc_t;
    std::vector<Timer> norm_t;
    OperatorStatistics<D> operStat;

    MWNodeVector<D> makeOperBand(const MWNode<D> &gNode, std::vector<NodeIndex<D>> &idx_band);

    void initTimers();
    void clearTimers();
    void printTimers() const;

    void calcNode(MWNode<D> &node) override;
    void postProcess() override {
        printTimers();
        clearTimers();
        initTimers();
    }

    void applyOperator(OperatorState<D> &os);
    void applyOperator_bw0(OperatorState<D> &os);
    void tensorApplyOperComp(OperatorState<D> &os);
};

} // namespace mrcpp
