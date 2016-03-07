/*
 *
 *
 *  \date Jul 8, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */

#include "FilterCache.h"

using namespace std;
using namespace Eigen;

template <int T>
FilterCache<T>::FilterCache() {
    switch (T) {
    case (Interpol):
        type = Interpol;
        break;
    case (Legendre):
        type = Legendre;
        break;
    default:
        MSG_ERROR("Invalid filter type: " << T);
    }
    this->libPath = MWFilter::getDefaultLibrary();
}

template <int T>
void FilterCache<T>::load(int order) {
    SET_CACHE_LOCK();
    if (not hasId(order)) {
        MWFilter *f = new MWFilter(order, type, this->libPath);
        int memo = f->getFilter().size() * sizeof(double);
        ObjectCache<MWFilter>::load(order, f, memo);
    }
    UNSET_CACHE_LOCK();
}

template <int T>
MWFilter &FilterCache<T>::get(int order) {
    if (not hasId(order)) {
        load(order);
    }
    return ObjectCache<MWFilter>::get(order);
}

template <int T>
const MatrixXd &FilterCache<T>::getFilterMatrix(int order) {
    if (not hasId(order)) {
        load(order);
    }
    return ObjectCache<MWFilter>::get(order).getFilter();
}

template class FilterCache<Interpol>;
template class FilterCache<Legendre>;
