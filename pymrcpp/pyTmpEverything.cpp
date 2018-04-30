/*
 *  \date Apr 26, 2018
 *  \author Magnar Bjørgve <magnar.bjorgve@uit.no> \n
 *          Hylleraas Centre for Quantum Molecular Sciences \n
 *          UiT - The Arctic University of Norway
 */

#include "PyBoundingBox.h"

#include "trees/BoundingBox.h"
#include "trees/MultiResolutionAnalysis.h"
#include "trees/MWTree.h"
#include "trees/FunctionTree.h"
#include "trees/FunctionTreeVector.h"

#include "core/ScalingBasis.h"
#include "core/InterpolatingBasis.h"
#include "core/LegendreBasis.h"

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"
#include "pybind11/functional.h"
#include "pybind11/eigen.h"

using namespace mrcpp;
namespace py = pybind11;


void pyTmpEverything(py::module &m) {

    std::stringstream funcTreeVecName;
    funcTreeVecName << "FunctionTreeVector" << 3 << "D";
    py::class_<FunctionTreeVector<3>> (m, funcTreeVecName.str().data())
            .def(py::init<>())
            .def("size", &FunctionTreeVector<3>::size)
            .def("push_back", py::overload_cast<double, FunctionTree<3> *>(&FunctionTreeVector<3>::push_back))
            .def("push_back", py::overload_cast<FunctionTree<3> *>(&FunctionTreeVector<3>::push_back));

    std::stringstream pyBoundBoxName;
    pyBoundBoxName << "BoundingBox" << 3 << "D";
    py::class_<PyBoundingBox<3>> (m, pyBoundBoxName.str().data())
            .def(py::init<int, py::array_t<int>, py::array_t <int>>())
            .def(py::init<int, int *,  int *>()) //1D cases can be initialized without array type input
            .def("getScale", &PyBoundingBox<3>::getScale);

    std::stringstream multResAnaName;
    multResAnaName << "MultiResolutionAnalysis" << 3 << "D";
    py::class_<MultiResolutionAnalysis<3>> (m, multResAnaName.str().data())
            .def(py::init<PyBoundingBox<3>, ScalingBasis, int>()) // Seems to work with PyBoundingBox here
            .def("getOrder", &MultiResolutionAnalysis<3>::getOrder)
            .def("getMaxDepth", &MultiResolutionAnalysis<3>::getMaxDepth)
            .def("getMaxScale", &MultiResolutionAnalysis<3>::getMaxScale);

    std::stringstream mwTreeName;
    mwTreeName << "MWTree" << 3 << "D";
    py::class_<MWTree<3>> mwtree(m, mwTreeName.str().data());
    mwtree
            .def(py::init<MultiResolutionAnalysis<3>>())
            .def("getSquareNorm", &MWTree<3>::getSquareNorm);

    std::stringstream funcTreeName;
    funcTreeName << "FunctionTree" << 3 << "D";
    py::class_<FunctionTree<3>> (m, funcTreeName.str().data(), mwtree)
            .def(py::init<MultiResolutionAnalysis<3>>())
            .def("integrate", &FunctionTree<3>::integrate)
            .def("clear", &FunctionTree<3>::clear)
            .def("normalize", &FunctionTree<3>::normalize);
            //.def("evalf", py::overload_cast<double>(&FunctionTree<3>::evalf))
            //.def("evalf", py::overload_cast<double, double>(&FunctionTree<3>::evalf))
            //.def("evalf", py::overload_cast<double, double, double>(&FunctionTree<3>::evalf));


    py::class_<ScalingBasis> scalingbasis(m, "ScalingBasis");
         scalingbasis.def(py::init<int, int>());

    py::class_<InterpolatingBasis> (m, "InterpolatingBasis", scalingbasis)
        .def(py::init<int>())
        .def("getScalingOrder", &InterpolatingBasis::getScalingOrder);

    py::class_<LegendreBasis> (m, "LegendreBasis", scalingbasis)
        .def(py::init<int>())
        .def("getScalingOrder", &LegendreBasis::getScalingOrder);


}