#ifndef WAVELETADAPTOR_H
#define WAVELETADAPTOR_H

#include "TreeAdaptor.h"
#include "constants.h"

template<int D>
class WaveletAdaptor : public TreeAdaptor<D> {
public:
    WaveletAdaptor(double pr,
		   int ms,
		   bool ap = false,
		   double sf = 1.0)
            : TreeAdaptor<D>(ms),
              prec(pr),
              absPrec(ap),
	      splitFac(sf) { }
    virtual ~WaveletAdaptor() { }

    void setPrecision(double pr) { this->prec = pr; }
    void setAbsPrec(bool ap) { this->absPrec = ap; }

protected:
    double prec;
    bool absPrec;
    double splitFac;

    virtual bool splitNode(const MWNode<D> &node) const;
    double getWaveletThreshold(double norm, int scale) const;
};

#endif // WAVELETADAPTOR_H
