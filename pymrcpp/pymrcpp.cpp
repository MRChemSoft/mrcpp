/*
 *  \date Mar 08, 2018
 *  \author Magnar Bjørgve <magnar.bjorgve@uit.no> \n
 *          Hylleraas Centre for Quantum Molecular Sciences \n
 *          UiT - The Arctic University of Norway
 */


#include "pybind11/pybind11.h"

namespace py = pybind11;

//template <int D> void pyFundamentalTypes(py::module &);
//template <int D> void pyGauss(py::module &);
void pyBases(py::module &);
template<int D> void pyMethods(py::module &);
//void pyProject1D(py::module &);
//void pyProject2D(py::module &);
void pyProject3D(py::module &);

template<int D> void pyOperators(py::module &);

void pyTmpEverything(py::module &);

PYBIND11_MODULE(pymrcpp, m) {

  //  pyProject1D(m);
  //  pyProject2D(m);
    pyProject3D(m);

    pyTmpEverything(m);

    //pyFundamentalTypes<1>(m);
    //pyFundamentalTypes<2>(m);
    //pyFundamentalTypes<3>(m);

    pyBases(m);
//
//    pyMethods<1>(m);
//    pyMethods<2>(m);
    pyMethods<3>(m);
//
//    pyOperators<1>(m);
//    pyOperators<2>(m);
    pyOperators<3>(m);
//
    //pyGauss<1>(m);
    //pyGauss<2>(m);
    //pyGauss<3>(m);

}
