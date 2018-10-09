/*
 *  \date Apr 25, 2018
 *  \author Magnar Bjørgve <magnar.bjorgve@uit.no> \n
 *          Hylleraas Centre for Quantum Molecular Sciences \n
 *          UiT - The Arctic University of Norway
 */


#include "PyBoundingBox.h"
using namespace mrcpp;

namespace py = pybind11;

template <int D>
PyBoundingBox<D>::PyBoundingBox(int n, py::array_t<int> l, py::array_t<int> nb)
    : BoundingBox<D>() {

    auto bufl = l.request();
    auto bufnb = nb.request();

    const int *lPtr = (const int *) bufl.ptr;
    const int *nbPtr = (const int *) bufnb.ptr;

    this->cornerIndex = NodeIndex<D>(n, lPtr);
    this->setNBoxes(nbPtr);
    this->setDerivedParameters();
}

template class PyBoundingBox<1>;
template class PyBoundingBox<2>;
template class PyBoundingBox<3>;
