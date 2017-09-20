/*
 *
 *
 *  \date Jul 18, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */

#include "CrossCorrelationCache.h"
#include "constants.h"

using namespace std;
using namespace Eigen;

template <int T>
CrossCorrelationCache<T>::CrossCorrelationCache() {
    switch (T) {
    case (Interpol):
        type = Interpol;
        break;
    case (Legendre):
        type = Legendre;
        break;
    default:
        MSG_ERROR("Invalid CrossCorrelation type: " << T)
    }
    this->libPath = CrossCorrelation::getDefaultLibrary();
}

template <int T>
void CrossCorrelationCache<T>::load(int order) {
    SET_CACHE_LOCK();
    if (not hasId(order)) {
        CrossCorrelation *ccc = new CrossCorrelation(order, type, this->libPath);
        int memo = ccc->getLMatrix().size() * 2 * sizeof(double);
        ObjectCache<CrossCorrelation>::load(order, ccc, memo);
    }
    UNSET_CACHE_LOCK();
}

template <int T>
CrossCorrelation &CrossCorrelationCache<T>::get(int order) {
    if (not hasId(order)) {
        load(order);
    }
    return ObjectCache<CrossCorrelation>::get(order);
}

template <int T>
const Eigen::MatrixXd &CrossCorrelationCache<T>::getLMatrix(int order) {
    if (not hasId(order)) {
        load(order);
    }
    return ObjectCache<CrossCorrelation>::get(order).getLMatrix();
}

template <int T>
const Eigen::MatrixXd &CrossCorrelationCache<T>::getRMatrix(int order) {
    if (not hasId(order)) {
        load(order);
    }
    return ObjectCache<CrossCorrelation>::get(order).getRMatrix();
}

template class CrossCorrelationCache<Interpol>;
template class CrossCorrelationCache<Legendre>;

