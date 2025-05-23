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
#include "trees/FunctionTreeVector.h"

#include "MRCPP/mrcpp_declarations.h"

namespace mrcpp {

template <int D, typename T> class ConvolutionCalculator final : public TreeCalculator<D, T> {
public:
    ConvolutionCalculator(double p, ConvolutionOperator<D> &o, FunctionTree<D, T> &f, int depth = MaxDepth);
    ~ConvolutionCalculator() override;

    MWNodeVector<D, T> *getInitialWorkVector(MWTree<D, T> &tree) const override;

    void setPrecFunction(const std::function<double(const NodeIndex<D> &idx)> &prec_func) { this->precFunc = prec_func; }
    void startManipulateOperator(bool excUnit) {
        this->manipulateOperator = true;
        this->onUnitcell = excUnit;
    }

private:
    int maxDepth;
    double prec;
    bool manipulateOperator{false};
    bool onUnitcell{false};
    ConvolutionOperator<D> *oper;
    FunctionTree<D, T> *fTree;
    std::vector<Timer *> band_t;
    std::vector<Timer *> calc_t;
    std::vector<Timer *> norm_t;

    OperatorStatistics operStat;
    std::vector<Eigen::MatrixXi *> bandSizes;
    std::function<double(const NodeIndex<D> &idx)> precFunc = [](const NodeIndex<D> &idx) { return 1.0; };

    static const int nComp = (1 << D);
    static const int nComp2 = (1 << D) * (1 << D);

    MWNodeVector<D, T> *makeOperBand(const MWNode<D, T> &gNode, std::vector<NodeIndex<D>> &idx_band);
    void fillOperBand(MWNodeVector<D, T> *band, std::vector<NodeIndex<D>> &idx_band, NodeIndex<D> &idx, const int *nbox, int dim);

    void initTimers();
    void clearTimers();
    void printTimers() const;

    void initBandSizes();
    int getBandSizeFactor(int i, int depth, const OperatorState<D, T> &os) const {
        int k = os.gt * this->nComp + os.ft;
        return (*this->bandSizes[i])(depth, k);
    }

    void calcBandSizeFactor(Eigen::MatrixXi &bs, int depth, const BandWidth &bw);

    void calcNode(MWNode<D, T> &node) override;
    void postProcess() override {
        printTimers();
        clearTimers();
        initTimers();
    }

    void applyOperComp(OperatorState<D, T> &os);
    void applyOperator(int i, OperatorState<D, T> &os);
    void tensorApplyOperComp(OperatorState<D, T> &os);

    void touchParentNodes(MWTree<D, T> &tree) const;
};

} // namespace mrcpp
