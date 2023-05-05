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

/** OperatorState is a simple helper class for operator application.
 * It keeps track of various state dependent variables and memory
 * regions. We cannot have some of this information directly in OperatorFunc
 * because of multi-threading issues.
 */

#pragma once

#include <vector>

#include <Eigen/Core>

#include "trees/MWNode.h"
#include "utils/math_utils.h"

namespace mrcpp {

#define GET_OP_IDX(FT, GT, ID) (2 * ((GT >> ID) & 1) + ((FT >> ID) & 1))

template <int D> class OperatorState final {
public:
    OperatorState(MWNode<D> &gn, double *scr1)
            : gNode(&gn) {
        this->kp1 = this->gNode->getKp1();
        this->kp1_d = this->gNode->getKp1_d();
        this->kp1_2 = math_utils::ipow(this->kp1, 2);
        this->kp1_dm1 = math_utils::ipow(this->kp1, D - 1);
        this->gData = this->gNode->getCoefs();
        this->maxDeltaL = -1;

        double *scr2 = scr1 + this->kp1_d;

        for (int i = 1; i < D; i++) {
            if (IS_ODD(i)) {
                this->aux[i] = scr2;
            } else {
                this->aux[i] = scr1;
            }
        }
    }

    OperatorState(MWNode<D> &gn, std::vector<double> scr1)
            : OperatorState(gn, scr1.data()) {}
    void setFNode(MWNode<D> &fn) {
        this->fNode = &fn;
        this->fData = this->fNode->getCoefs();
    }
    void setFIndex(NodeIndex<D> &idx) {
        this->fIdx = &idx;
        calcMaxDeltaL();
    }
    void setGComponent(int gt) {
        this->aux[D] = this->gData + gt * this->kp1_d;
        this->gt = gt;
    }
    void setFComponent(int ft) {
        this->aux[0] = this->fData + ft * this->kp1_d;
        this->ft = ft;
    }

    int getMaxDeltaL() const { return this->maxDeltaL; }
    int getOperIndex(int i) const { return GET_OP_IDX(this->ft, this->gt, i); }

    double **getAuxData() { return this->aux; }
    double **getOperData() { return this->oData; }

    friend class ConvolutionCalculator<D>;
    friend class DerivativeCalculator<D>;

private:
    int ft;
    int gt;
    int maxDeltaL;
    double fThreshold;
    double gThreshold;
    // Shorthands
    int kp1;
    int kp1_2;
    int kp1_d;
    int kp1_dm1;

    MWNode<D> *gNode;
    MWNode<D> *fNode;
    NodeIndex<D> *fIdx;

    double *aux[D + 1];
    double *gData;
    double *fData;
    double *oData[D];

    void calcMaxDeltaL() {
        const auto &gl = this->gNode->getNodeIndex();
        const auto &fl = *this->fIdx;
        int max_dl = 0;
        for (int d = 0; d < D; d++) {
            int dl = std::abs(fl[d] - gl[d]);
            if (dl > max_dl) { max_dl = dl; }
        }
        this->maxDeltaL = max_dl;
    }
};

} // namespace mrcpp
