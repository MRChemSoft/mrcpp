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

/*
 *
 *
 *  \date Jul 18, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */

#include "CrossCorrelationCache.h"
#include "utils/Printer.h"

#include "MRCPP/constants.h"

using namespace Eigen;

namespace mrcpp {

template <int T> CrossCorrelationCache<T>::CrossCorrelationCache() {
    switch (T) {
        case (Interpol):
            type = Interpol;
            break;
        case (Legendre):
            type = Legendre;
            break;
        default:
            MSG_ERROR("Invalid CrossCorrelation type: " << T)
    }
}

template <int T> void CrossCorrelationCache<T>::load(int order) {
    SET_CACHE_LOCK();
    if (not hasId(order)) {
        auto *ccc = new CrossCorrelation(order, type);
        int memo = ccc->getLMatrix().size() * 2 * sizeof(double);
        ObjectCache<CrossCorrelation>::load(order, ccc, memo);
    }
    UNSET_CACHE_LOCK();
}

template <int T> CrossCorrelation &CrossCorrelationCache<T>::get(int order) {
    if (not hasId(order)) { load(order); }
    return ObjectCache<CrossCorrelation>::get(order);
}

template <int T> const Eigen::MatrixXd &CrossCorrelationCache<T>::getLMatrix(int order) {
    if (not hasId(order)) { load(order); }
    return ObjectCache<CrossCorrelation>::get(order).getLMatrix();
}

template <int T> const Eigen::MatrixXd &CrossCorrelationCache<T>::getRMatrix(int order) {
    if (not hasId(order)) { load(order); }
    return ObjectCache<CrossCorrelation>::get(order).getRMatrix();
}

template class CrossCorrelationCache<Interpol>;
template class CrossCorrelationCache<Legendre>;

} // namespace mrcpp
