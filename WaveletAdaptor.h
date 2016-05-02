#ifndef WAVELETADAPTOR_H
#define WAVELETADAPTOR_H

#include "TreeAdaptor.h"
#include "constants.h"

template<int D>
class WaveletAdaptor : public TreeAdaptor<D> {
public:
    WaveletAdaptor(double p = -1.0, int n = MaxScale, bool a = false)
            : prec(p),
              absPrec(a),
              maxScale(n) { }
    WaveletAdaptor(const WaveletAdaptor<D> &adap)
            : prec(adap.prec),
              absPrec(adap.absPrec),
              maxScale(adap.maxScale) { }
    virtual ~WaveletAdaptor() { }
    virtual TreeAdaptor<D> *copy() const { return new WaveletAdaptor<D>(*this); }

protected:
    const double prec;
    const bool absPrec;
    const int maxScale;

    virtual bool splitNode(const MWNode<D> &node) const;
    double getWaveletThreshold(double norm, int scale) const;
};

#endif // WAVELETADAPTOR_H
