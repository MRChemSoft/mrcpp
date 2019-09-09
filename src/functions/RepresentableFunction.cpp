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

/**
 *
 *
 *  \date April 30, 2010
 *  \author Stig Rune Jensen  \n
 *          CTCC, University of Tromsø
 *
 *
 */

#include "RepresentableFunction.h"
#include "utils/Printer.h"

namespace mrcpp {

template <int D> RepresentableFunction<D>::RepresentableFunction(const double *a, const double *b) {
    if (a == nullptr or b == nullptr) {
        this->bounded = false;
        this->A = nullptr;
        this->B = nullptr;
    } else {
        this->bounded = true;
        this->A = new double[D];
        this->B = new double[D];
        for (int d = 0; d < D; d++) {
            if (a[d] > b[d]) { MSG_ERROR("Lower bound > Upper bound."); }
            this->A[d] = a[d];
            this->B[d] = b[d];
        }
    }
}

/** Constructs a new function with same bounds as the input function */
template <int D> RepresentableFunction<D>::RepresentableFunction(const RepresentableFunction<D> &func) {
    if (func.isBounded()) {
        this->bounded = true;
        this->A = new double[D];
        this->B = new double[D];
        for (int d = 0; d < D; d++) {
            A[d] = func.getLowerBounds()[d];
            B[d] = func.getUpperBounds()[d];
        }
    } else {
        this->bounded = false;
        this->A = nullptr;
        this->B = nullptr;
    }
}

/** Copies function, not bounds. Use copy constructor if you want an
 * identical function. */
template <int D> RepresentableFunction<D> &RepresentableFunction<D>::operator=(const RepresentableFunction<D> &func) {
    return *this;
}

template <int D> RepresentableFunction<D>::~RepresentableFunction() {
    if (this->isBounded()) {
        delete[] this->A;
        delete[] this->B;
    }
    this->A = nullptr;
    this->B = nullptr;
}

template <int D> void RepresentableFunction<D>::setBounds(const double *a, const double *b) {
    if (a == nullptr or b == nullptr) { MSG_ERROR("Invalid arguments"); }
    if (not isBounded()) {
        this->bounded = true;
        this->A = new double[D];
        this->B = new double[D];
    }
    for (int d = 0; d < D; d++) {
        if (a[d] > b[d]) { MSG_ERROR("Lower bound > Upper bound."); }
        this->A[d] = a[d];
        this->B[d] = b[d];
    }
}

template <int D> bool RepresentableFunction<D>::outOfBounds(const Coord<D> &r) const {
    if (not isBounded()) { return false; }
    for (int d = 0; d < D; d++) {
        if (r[d] < getLowerBound(d)) { return true; }
        if (r[d] > getUpperBound(d)) { return true; }
    }
    return false;
}

template class RepresentableFunction<1>;
template class RepresentableFunction<2>;
template class RepresentableFunction<3>;

} // namespace mrcpp
