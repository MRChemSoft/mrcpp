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

#include "TreeCalculator.h"
#include "operators/OperatorStatistics.h"

#include "mrcpp_declarations.h"

namespace mrcpp {

template<int D>
class ConvolutionCalculator final : public TreeCalculator<D> {
public:
    ConvolutionCalculator(double p, ConvolutionOperator<D> &o, FunctionTree<D> &f, int depth = MaxDepth);
    ~ConvolutionCalculator();

    MWNodeVector<D>* getInitialWorkVector(MWTree<D> &tree) const;

private:
    int maxDepth;
    double prec;
    ConvolutionOperator<D> *oper;
    FunctionTree<D> *fTree;
    std::vector<Timer *> band_t;
    std::vector<Timer *> calc_t;
    std::vector<Timer *> norm_t;

    OperatorStatistics<D> operStat;
    std::vector<Eigen::MatrixXi *> bandSizes;

    static const int nComp = (1 << D);
    static const int nComp2 = (1 << D) * (1 << D);

    MWNodeVector<D>* makeOperBand(const MWNode<D> &gNode, std::vector<NodeIndex<D> > &idx_band);
    void fillOperBand(MWNodeVector<D> *band, std::vector<NodeIndex<D> > &idx_band, NodeIndex<D> &idx, const int *nbox, int dim);

    void initTimers();
    void clearTimers();
    void printTimers() const;

    void initBandSizes();
    int getBandSizeFactor(int i, int depth,const OperatorState<D> &os) const;
    void calcBandSizeFactor(Eigen::MatrixXi &bs, int depth, const BandWidth &bw);

    void calcNode(MWNode<D> &node);
    void postProcess() {
        printTimers();
        clearTimers();
        initTimers();
    }

    void applyOperComp(OperatorState<D> &os);
    void applyOperator(OperatorState<D> &os);
    void tensorApplyOperComp(OperatorState<D> &os);
};

}
