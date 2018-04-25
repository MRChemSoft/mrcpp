/*
 *  \date Mar 22, 2018
 *  \author Magnar Bjørgve <magnar.bjorgve@uit.no> \n
 *          Hylleraas Centre for Quantum Molecular Sciences \n
 *          UiT - The Arctic University of Norway
 */
//#include "pybind11/pybind11.h"
//#include "pybind11/numpy.h"

#include "mwfunctions/RepresentableFunction.h"
#include "mwfunctions/Gaussian.h"
#include "mwfunctions/GaussFunc.h"
#include "PyRepresentableFunction.h" // Trampoline Class for Representable Function

using namespace mrcpp;
namespace py = pybind11;

template<int D>
void pyGauss(py::module &m) {

    std::stringstream repFuncName;
    repFuncName << "RepresentableFunction" << D << "D";
    py::class_<RepresentableFunction<D>, PyRepresentableFunction<D>> repfunc(m, repFuncName.str().data());
    repfunc
        .def(py::init<>())
        .def("evalf", &RepresentableFunction<D>::evalf);

    std::stringstream gaussianName;
    gaussianName << "Gaussian" << D << "D";
    py::class_<Gaussian<D>> gaussian(m, gaussianName.str().data(), repfunc);

    std::stringstream gausFuncName;
    gausFuncName << "GaussFunc" << D << "D";
    py::class_<GaussFunc<D>>(m, gausFuncName.str().data(), gaussian)
        .def(py::init<double, double, py::array_t <double>, py::array_t <double>>())
        .def("evalf", py::overload_cast<py::array_t <double>>(&GaussFunc<D>::evalf))
        .def("calcCoulombEnergy", &GaussFunc<D>::calcCoulombEnergy);
}

template void pyGauss<1>(py::module &m);
template void pyGauss<2>(py::module &m);
template void pyGauss<3>(py::module &m);
