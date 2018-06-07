/*
 *  \date Apr 26, 2018
 *  \author Magnar Bjørgve <magnar.bjorgve@uit.no> \n
 *          Hylleraas Centre for Quantum Molecular Sciences \n
 *          UiT - The Arctic University of Norway
 */

#include "treebuilders/project.h"
#include "trees/FunctionTree.h"
#include "PyAnalyticFunction.h"

#include "pybind11/functional.h"
#include "pybind11/pybind11.h"
#include "pyProject.h"
#include <utility>

namespace py = pybind11;
using namespace mrcpp;

void project3D(double prec,
               FunctionTree<3> &out,
               std::function<double (double x, double y, double z)> func,
               int maxIter) {
    PyAnalyticFunction3D inp(std::move(func));
    project(prec, out, inp, maxIter);
}

void project2D(double prec,
                      FunctionTree<2> &out,
                      std::function<double (double x, double y)> func,
                      int maxIter) {
    PyAnalyticFunction2D inp(std::move(func));
    project(prec, out, inp, maxIter);
}

void project1D(double prec,
                      FunctionTree<1> &out,
                      std::function<double (double x)> func,
                      int maxIter) {
    PyAnalyticFunction1D inp(std::move(func));
    project(prec, out, inp, maxIter);
}

void pyProject1D(py::module &m) {
    m.def("project", py::overload_cast<double, FunctionTree<1> &, std::function<double (double)>, int>(&project1D),
          py::arg("prec"), py::arg("out"), py::arg("func"), py::arg("maxIter")= -1);
}

void pyProject2D(py::module &m) {
    m.def("project", py::overload_cast<double, FunctionTree<2> &, std::function<double (double, double)>, int>(&project2D),
          py::arg("prec"), py::arg("out"), py::arg("func"), py::arg("maxIter")= -1);
}

void pyProject3D(py::module &m) {
    m.def("project", py::overload_cast<double, FunctionTree<3> &, std::function<double (double, double, double)>, int>(&project3D),
          py::arg("prec"), py::arg("out"), py::arg("func"), py::arg("maxIter")= -1);
}
