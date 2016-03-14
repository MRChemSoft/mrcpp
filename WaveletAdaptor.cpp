#include "WaveletAdaptor.h"
#include "MWTree.h"
#include "MWNode.h"

#include "constants.h"

using namespace std;

template<int D>
bool WaveletAdaptor<D>::splitNode(const MWNode<D> &node) const {
    if (this->prec < 0.0) {
        return false;
    }
    int scale = node.getScale();
    if ((scale + 1) >= this->maxScale) {
        println(10, "Maximum scale reached: " << scale + 1);
        return false;
    }
    double t_norm = 1.0;
    double sq_norm = node.getMWTree().getSquareNorm();
    if (not (sq_norm > 0.0) and not this->absPrec) {
        t_norm = sqrt(sq_norm);
    }
    double w_norm = node.getWaveletNorm();
    double thrs = getWaveletThreshold(t_norm, scale);

    if (w_norm > thrs) {
        return true;
    } else {
        return false;
    }
}

/** Calculate the threshold for the wavelet norm.
  *
  * Calculates the threshold that has to be met in the wavelet norm in order to
  * guarantee the precision in the function representation. Depends on the
  * square norm of the function and the requested relative accuracy. */
template<int D>
double WaveletAdaptor<D>::getWaveletThreshold(double norm, int scale) const {
    double expo = (0.5 * (scale + 1));
    double thrs_1 = 2.0 * MachinePrec;
    double thrs_2 = norm * this->prec * pow(2.0, -expo);
    return max(thrs_1, thrs_2);
}

template class WaveletAdaptor<1>;
template class WaveletAdaptor<2>;
template class WaveletAdaptor<3>;
