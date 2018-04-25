/*
 *  \date Mar 16, 2018
 *  \author Magnar Bjørgve <magnar.bjorgve@uit.no> \n
 *          Hylleraas Centre for Quantum Molecular Sciences \n
 *          UiT - The Arctic University of Norway
 */
#pragma once
//#include "pybind11/pybind11.h"
#include "mwfunctions/RepresentableFunction.h"


template <int D>
class PyRepresentableFunction : public mrcpp::RepresentableFunction<D> {
public:
    using mrcpp::RepresentableFunction<D>::RepresentableFunction;

    double evalf(const double *r) const override {
        PYBIND11_OVERLOAD_PURE(
                double,
                mrcpp::RepresentableFunction<D>,
                evalf,
                r
        );
    }

};
