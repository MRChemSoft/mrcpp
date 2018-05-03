/*
 *  \date May 02, 2018
 *  \author Magnar Bjørgve <magnar.bjorgve@uit.no> \n
 *          Hylleraas Centre for Quantum Molecular Sciences \n
 *          UiT - The Arctic University of Norway
 */

#pragma once

#include "trees/FunctionTree.h"

template <int D>
class PyFunctionTree: public mrcpp::FunctionTree<3> {
public:
    using mrcpp::FunctionTree<3>::FunctionTree;
    virtual double evalf3D(double x, double y, double z);
};