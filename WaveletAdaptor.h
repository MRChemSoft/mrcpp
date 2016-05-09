#ifndef WAVELETADAPTOR_H
#define WAVELETADAPTOR_H

#include "TreeAdaptor.h"
#include "constants.h"

template<int D>
class WaveletAdaptor : public TreeAdaptor<D> {
public:
    WaveletAdaptor(double pr = -1.0, int ms = MaxScale, bool ap = false)
            : TreeAdaptor<D>(ms),
              prec(pr),
              absPrec(ap) { }
    virtual ~WaveletAdaptor() { }

    void setPrecision(double pr) { this->prec = pr; }
    void setAbsPrec(bool ap) { this->absPrec = ap; }

protected:
    double prec;
    bool absPrec;

    virtual bool splitNode(const MWNode<D> &node) const;
    double getWaveletThreshold(double norm, int scale) const;
};

#endif // WAVELETADAPTOR_H
