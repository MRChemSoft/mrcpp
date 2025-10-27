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

#include "OperatorStatistics.h"
#include "trees/MWNode.h"

using namespace Eigen;

namespace mrcpp {

OperatorStatistics::OperatorStatistics()
        : nThreads(mrcpp_get_max_threads())
        , totFCount(0)
        , totGCount(0)
        , totGenCount(0)
        , fCount(nullptr)
        , gCount(nullptr)
        , genCount(nullptr)
        , totCompCount(nullptr)
        , compCount(nullptr) {

    this->totCompCount = new Matrix<int, 8, 8>;
    this->totCompCount->setZero();

    this->fCount = new int[this->nThreads];
    this->gCount = new int[this->nThreads];
    this->genCount = new int[this->nThreads];
    this->compCount = new Matrix<int, 8, 8> *[this->nThreads];
    for (int i = 0; i < this->nThreads; i++) {
        this->compCount[i] = new Matrix<int, 8, 8>;
        this->compCount[i]->setZero();
        this->fCount[i] = 0;
        this->gCount[i] = 0;
        this->genCount[i] = 0;
    }
}

OperatorStatistics::~OperatorStatistics() {
    for (int i = 0; i < this->nThreads; i++) { delete this->compCount[i]; }
    delete[] this->compCount;
    delete[] this->fCount;
    delete[] this->gCount;
    delete[] this->genCount;
    delete totCompCount;
}

/** Sum all node counters from all threads. */
void OperatorStatistics::flushNodeCounters() {
    for (int i = 0; i < this->nThreads; i++) {
        this->totFCount += this->fCount[i];
        this->totGCount += this->gCount[i];
        this->totGenCount += this->genCount[i];
        *this->totCompCount += *this->compCount[i];
        this->fCount[i] = 0;
        this->gCount[i] = 0;
        this->genCount[i] = 0;
        this->compCount[i]->setZero();
    }
}

/** Increment g-node usage counter. Needed for load balancing. */
template <int D, typename T> void OperatorStatistics::incrementGNodeCounters(const MWNode<D, T> &gNode) {
    int thread = mrcpp_get_thread_num();
    this->gCount[thread]++;
}

/** Increment operator application counter. */
template <int D, typename T> void OperatorStatistics::incrementFNodeCounters(const MWNode<D, T> &fNode, int ft, int gt) {
    int thread = mrcpp_get_thread_num();
    this->fCount[thread]++;
    (*this->compCount[thread])(ft, gt) += 1;
    if (fNode.isGenNode()) { this->genCount[thread]++; }
}

std::ostream &OperatorStatistics::print(std::ostream &o) const {
    o << std::setw(8);
    o << "*OperatorFunc statistics: " << std::endl << std::endl;
    o << "  Total calculated gNodes      : " << this->totGCount << std::endl;
    o << "  Total applied fNodes         : " << this->totFCount << std::endl;
    o << "  Total applied genNodes       : " << this->totGenCount << std::endl << std::endl;
    o << "  By components:" << std::endl << *this->totCompCount << std::endl;
    return o;
}

template void OperatorStatistics::incrementFNodeCounters<1, double>(const MWNode<1, double> &fNode, int ft, int gt);
template void OperatorStatistics::incrementFNodeCounters<2, double>(const MWNode<2, double> &fNode, int ft, int gt);
template void OperatorStatistics::incrementFNodeCounters<3, double>(const MWNode<3, double> &fNode, int ft, int gt);
template void OperatorStatistics::incrementFNodeCounters<1, ComplexDouble>(const MWNode<1, ComplexDouble> &fNode, int ft, int gt);
template void OperatorStatistics::incrementFNodeCounters<2, ComplexDouble>(const MWNode<2, ComplexDouble> &fNode, int ft, int gt);
template void OperatorStatistics::incrementFNodeCounters<3, ComplexDouble>(const MWNode<3, ComplexDouble> &fNode, int ft, int gt);
template void OperatorStatistics::incrementGNodeCounters<1, double>(const MWNode<1, double> &gNode);
template void OperatorStatistics::incrementGNodeCounters<2, double>(const MWNode<2, double> &gNode);
template void OperatorStatistics::incrementGNodeCounters<3, double>(const MWNode<3, double> &gNode);
template void OperatorStatistics::incrementGNodeCounters<1, ComplexDouble>(const MWNode<1, ComplexDouble> &gNode);
template void OperatorStatistics::incrementGNodeCounters<2, ComplexDouble>(const MWNode<2, ComplexDouble> &gNode);
template void OperatorStatistics::incrementGNodeCounters<3, ComplexDouble>(const MWNode<3, ComplexDouble> &gNode);

} // namespace mrcpp
