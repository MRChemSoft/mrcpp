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
    double t_norm = 1.0;
    double sq_norm = node.getMWTree().getSquareNorm();
    if (not (sq_norm > 0.0) and not this->absPrec) {
        t_norm = sqrt(sq_norm);
    }
    int scale = node.getScale();
    double thrs = getWaveletThreshold(t_norm, scale);

    double w_norm = sqrt(node.getWaveletNorm());
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
    double scale_fac = 1.0;
    if (this->splitFac > MachineZero) {
    	double expo = 0.5 * this->splitFac * (scale + 1);
	scale_fac = pow(2.0, -expo);
    }
    double thrs_1 = 2.0 * MachinePrec;
    double thrs_2 = norm * this->prec * scale_fac;
    return max(thrs_1, thrs_2);
}

template class WaveletAdaptor<1>;
template class WaveletAdaptor<2>;
template class WaveletAdaptor<3>;
