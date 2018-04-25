/*
 *  \date Mar 20, 2018
 *  \author Magnar Bjørgve <magnar.bjorgve@uit.no> \n
 *          Hylleraas Centre for Quantum Molecular Sciences \n
 *          UiT - The Arctic University of Norway
 */

#include "pybind11/pybind11.h"
//#include "pybind11/functional.h"

#include "mwtrees/FunctionTree.h"
#include "mwtrees/FunctionTreeVector.h"
#include "mwoperators/DerivativeOperator.h"
#include "mwoperators/ConvolutionOperator.h"
#include "mwbuilders/add.h"
#include "mwbuilders/multiply.h"
#include "mwbuilders/apply.h"
#include "mwbuilders/grid.h"
#include "mwbuilders/project.h"

namespace py = pybind11;
using namespace mrcpp;
/*
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
*/
template <int D>
void pyMethods(py::module &m) {


    m.def("add", py::overload_cast<double, FunctionTree<D> &, double , FunctionTree<D> &, double, FunctionTree<D> &, int>(&add<D>),
        py::arg("prec"), py::arg("out"), py::arg("a"),  py::arg("three_a"), py::arg("b"), py::arg("three_b"), py::arg("maxIter") = -1,
        "Adds to function threes");

    m.def("project", py::overload_cast<double, FunctionTree<D> &, RepresentableFunction<D> &, int>(&project<D>),
        py::arg("prec"), py::arg("out"), py::arg("inp"), py::arg("maxIter")= -1);

    m.def("multiply", py::overload_cast<double, FunctionTree<D> &, double, FunctionTree<D> &, FunctionTree<D> &, int >(&multiply<D>));
    m.def("multiply", py::overload_cast<double, FunctionTree<D> &, FunctionTreeVector<D> &, int >(&multiply<D>));

    m.def("build_grid", &build_grid<D>);
    m.def("copy_grid", py::overload_cast<FunctionTree<D> &, FunctionTree<D> &, int >(&copy_grid<D>));
    m.def("copy_grid", py::overload_cast<FunctionTree<D> &, FunctionTreeVector<D> &, int>(&copy_grid<D>));
    m.def("clear_grid", &clear_grid<D>);

    m.def("apply", py::overload_cast<double, FunctionTree<D> &, ConvolutionOperator<D> &, FunctionTree<D> &, int>(&apply<D>),
        py::arg("prec"), py::arg("out"), py::arg("oper"), py::arg("inp"), py::arg("maxIter") = -1);
    m.def("apply", py::overload_cast<FunctionTree<D> &, DerivativeOperator<D> &, FunctionTree<D> &, int>(&apply<D>));

    m.def("dot", &dot<D>);
}

template void pyMethods<1>(py::module &m);
template void pyMethods<2>(py::module &m);
template void pyMethods<3>(py::module &m);
