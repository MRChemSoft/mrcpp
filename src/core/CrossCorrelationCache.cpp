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
    MRCPP_SET_OMP_LOCK();
    if (not hasId(order)) {
        auto *ccc = new CrossCorrelation(order, type);
        int memo = ccc->getLMatrix().size() * 2 * sizeof(double);
        ObjectCache<CrossCorrelation>::load(order, ccc, memo);
    }
    MRCPP_UNSET_OMP_LOCK();
}

template <int T> CrossCorrelation &CrossCorrelationCache<T>::get(int order) {
    if (not hasId(order)) { load(order); }
    return ObjectCache<CrossCorrelation>::get(order);
}

template <int T> const Eigen::MatrixXd &CrossCorrelationCache<T>::getLMatrix(int order) {
    if (not hasId(order)) { load(order); }
    return ObjectCache<CrossCorrelation>::get(order).getLMatrix();
}


/** @brief Fetches the cross correlation coefficients.
 *
 * @param[in] order: Dimension of \f$ V_0 \subset L^2(\mathbb R) \f$ minus one,
 * that is the maximum degree @f( k @f) of polynomials in \f$ V_0 \subset L^2(0, 1) \f$.
 * @returns The right matrix of cross correlation coefficients.
 *
 * @details The cross correlation coefficients
 * \f[
 *   C^{(+)}_{ijp}
 *   =
 *   \int_0^1 dz
 *   \int_0^1 dx
 *   \phi_i(x)
 *   \phi_j(x - z)
 *   \phi_p(z)
 * \f]
 * with $i, j = 0, \ldots, k$ and $p = 0, \ldots, 2k + 1$.
 * They are grouped in the so called right matrix
 * \f[
 *   \begin{pmatrix}
 *       C^{(+)}_{000}
 *       &
 *       C^{(+)}_{001}
 *       &
 *       \ldots
 *       &
 *       C^{(+)}_{00,2k+1}
 *       \\
 *       C^{(+)}_{010}
 *       &
 *       C^{(+)}_{011}
 *       &
 *       \ldots
 *       &
 *       C^{(+)}_{01,2k+1}
 *       \\
 *       \ldots
 *       &
 *       \ldots
 *       &
 *       \ldots
 *       &
 *       \ldots
 *       \\
 *       C^{(+)}_{k, k - 1, 0}
 *       &
 *       C^{(+)}_{k, k - 1, 1}
 *       &
 *       \ldots
 *       &
 *       C^{(+)}_{k, k - 1, 2k+1}
 *       \\
 *       C^{(+)}_{kk0}
 *       &
 *       C^{(+)}_{kk1}
 *       &
 *       \ldots
 *       &
 *       C^{(+)}_{kk,2k+1}
 *   \end{pmatrix}
 * \f]
 * that is returned by the method.
 * 
 *
 */
template <int T> const Eigen::MatrixXd &CrossCorrelationCache<T>::getRMatrix(int order) {
    if (not hasId(order)) { load(order); }
    return ObjectCache<CrossCorrelation>::get(order).getRMatrix();
}

template class CrossCorrelationCache<Interpol>;
template class CrossCorrelationCache<Legendre>;

} // namespace mrcpp
