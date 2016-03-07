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

using namespace std;

template<int D>
RepresentableFunction<D>::RepresentableFunction(const double *a,
                                                const double *b) {
    if (a == 0 or b == 0) {
        this->bounded = false;
        this->A = 0;
        this->B = 0;
    } else {
        this->bounded = true;
        this->A = new double[D];
        this->B = new double[D];
        for (int d = 0; d < D; d++) {
            if (a[d] > b[d]) {
                MSG_ERROR("Lower bound > Upper bound.");
            }
            this->A[d] = a[d];
            this->B[d] = b[d];
        }
    }
}

/** Constructs a new function with same bounds as the input function */
template<int D>
RepresentableFunction<D>::RepresentableFunction(
        const RepresentableFunction<D> &func) {
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
        this->A = 0;
        this->B = 0;
    }
}

/** Copies function, not bounds. Use copy constructor if you want an
  * identical function. */
template<int D>
RepresentableFunction<D>& RepresentableFunction<D>::operator=(
        const RepresentableFunction<D> &func) {
    return *this;
}

template<int D>
RepresentableFunction<D>::~RepresentableFunction() {
    if (this->isBounded()) {
        delete[] this->A;
        delete[] this->B;
    }
    this->A = 0;
    this->B = 0;
}

template<int D>
void RepresentableFunction<D>::setBounds(const double *a, const double *b) {
    if (a == 0 or b == 0) {
        MSG_ERROR("Invalid arguments");
    }
    if (not isBounded()) {
        this->bounded = true;
        this->A = new double[D];
        this->B = new double[D];
    }
    for (int d = 0; d < D; d++) {
        if (a[d] > b[d]) {
            MSG_ERROR("Lower bound > Upper bound.");
        }
        this->A[d] = a[d];
        this->B[d] = b[d];
    }
}

template<int D>
bool RepresentableFunction<D>::outOfBounds(const double *r) const {
    if (not isBounded()) {
        return false;
    }
    for (int d = 0; d < D; d++) {
        if (r[d] < getLowerBound(d)) {
            return true;
        }
        if (r[d] > getUpperBound(d)) {
            return true;
        }
    }
    return false;
}

template class RepresentableFunction<1>;
template class RepresentableFunction<2>;
template class RepresentableFunction<3>;
