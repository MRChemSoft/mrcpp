/*
 *  \date May 02, 2018
 *  \author Magnar Bjørgve <magnar.bjorgve@uit.no> \n
 *          Hylleraas Centre for Quantum Molecular Sciences \n
 *          UiT - The Arctic University of Norway
 */

#include "PyFunctionTree.h"


template <int D>
double PyFunctionTree<D>::evalf3D(double x, double y, double z) {
    const double r[3] = {x, y, z};
    return mrcpp::FunctionTree<3>::evalf(r);
}

template class PyFunctionTree<3>;