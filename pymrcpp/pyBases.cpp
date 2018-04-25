/*
 *  \date Mar 08, 2018
 *  \author Magnar Bjørgve <magnar.bjorgve@uit.no> \n
 *          Hylleraas Centre for Quantum Molecular Sciences \n
 *          UiT - The Arctic University of Norway
 */

#include "pybind11/pybind11.h"

#include "mwcore/ScalingBasis.h"
#include "mwcore/InterpolatingBasis.h"
#include "mwcore/LegendreBasis.h"

using namespace mrcpp;
namespace py = pybind11;


void pyBases(py::module &m) {

py::class_<ScalingBasis> scalingbasis(m, "ScalingBasis");
    scalingbasis.def(py::init<int, int>());

py::class_<InterpolatingBasis> (m, "InterpolatingBasis", scalingbasis)
    .def(py::init<int>());

py::class_<LegendreBasis> (m, "LegendreBasis", scalingbasis)
    .def(py::init<int>());
}
