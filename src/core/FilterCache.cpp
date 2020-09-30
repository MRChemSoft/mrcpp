/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2020 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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
 *  \date Jul 8, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */

#include "FilterCache.h"
#include "utils/Printer.h"

#include "MRCPP/constants.h"

using namespace Eigen;

namespace mrcpp {

template <int T> FilterCache<T>::FilterCache() {
    switch (T) {
        case (Interpol):
            type = Interpol;
            break;
        case (Legendre):
            type = Legendre;
            break;
        default:
            MSG_ERROR("Invalid filter type: " << T);
    }
}

template <int T> void FilterCache<T>::load(int order) {
    MRCPP_SET_OMP_LOCK();
    if (not hasId(order)) {
        auto *f = new MWFilter(order, type);
        int memo = f->getFilter().size() * sizeof(double);
        ObjectCache<MWFilter>::load(order, f, memo);
    }
    MRCPP_UNSET_OMP_LOCK();
}

template <int T> MWFilter &FilterCache<T>::get(int order) {
    if (not hasId(order)) { load(order); }
    return ObjectCache<MWFilter>::get(order);
}

template <int T> const MatrixXd &FilterCache<T>::getFilterMatrix(int order) {
    if (not hasId(order)) { load(order); }
    return ObjectCache<MWFilter>::get(order).getFilter();
}

template class FilterCache<Interpol>;
template class FilterCache<Legendre>;

} // namespace mrcpp
