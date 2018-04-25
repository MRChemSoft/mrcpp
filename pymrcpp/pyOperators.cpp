/*
 *  \date Mar 21, 2018
 *  \author Magnar Bjørgve <magnar.bjorgve@uit.no> \n
 *          Hylleraas Centre for Quantum Molecular Sciences \n
 *          UiT - The Arctic University of Norway
 */


#include "pybind11/pybind11.h"

#include "mwtrees/MultiResolutionAnalysis.h"
#include "mwoperators/ABGVOperator.h"
#include "mwoperators/PoissonOperator.h"
#include "mwoperators/HelmholtzOperator.h"
#include "mwoperators/ConvolutionOperator.h"

using namespace mrcpp;
namespace py = pybind11;


template<int D>
void pyOperators(py::module &m) {

    std::stringstream DerOperatorName;
    DerOperatorName << "DerivativeOperator" << D << "D";
    py::class_<DerivativeOperator<D>> deriv(m, DerOperatorName.str().data());
    deriv
        .def(py::init<MultiResolutionAnalysis<D>>());

    std::stringstream ABGVOperatorName;
    ABGVOperatorName << "ABGVOperator" << D << "D";
    py::class_<ABGVOperator<D>> (m, ABGVOperatorName.str().data(), deriv)
        .def(py::init< MultiResolutionAnalysis<D> &, double, double >());


    std::stringstream ConvOperatorName;
    ConvOperatorName << "ConvolutionOperator" << D << "D";
    py::class_<ConvolutionOperator<D>> convop(m, ConvOperatorName.str().data());
    convop
        .def(py::init<MultiResolutionAnalysis<D> &, double>());

    if (D==3){
        py::class_<PoissonOperator> (m, "PoissonOperator", convop)
            .def(py::init<const MultiResolutionAnalysis<3> &, double >());
        py::class_<HelmholtzOperator> (m, "HelmholtzOperator", convop)
            .def(py::init<MultiResolutionAnalysis<3> &, double, double>());
    }

}

template void pyOperators<1>(py::module &m);
template void pyOperators<2>(py::module &m);
template void pyOperators<3>(py::module &m);

