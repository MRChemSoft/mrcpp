/*
 *  \date Apr 25, 2018
 *  \author Magnar Bjørgve <magnar.bjorgve@uit.no> \n
 *          Hylleraas Centre for Quantum Molecular Sciences \n
 *          UiT - The Arctic University of Norway
 */

#pragma once

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"

#include "mwtrees/BoundingBox.h"

template <int D>
class PyBoundingBox : public mrcpp::BoundingBox<D> {
public:
    using mrcpp::BoundingBox<D>::BoundingBox;
    PyBoundingBox(int, pybind11::array_t<int>, pybind11::array_t<int>);
};
