#pragma once

#include "TreeAdaptor.h"

namespace mrcpp {

template<int D>
class WaveletAdaptor : public TreeAdaptor<D> {
public:
    WaveletAdaptor(double pr, int ms, bool ap = false, double sf = 1.0)
        : TreeAdaptor<D>(ms),
          absPrec(ap),
          prec(pr),
          splitFac(sf) { }
    virtual ~WaveletAdaptor() = default;

protected:
    bool absPrec;
    double prec;
    double splitFac;

    virtual bool splitNode(const MWNode<D> &node) const {
        return node.splitCheck(this->prec, this->splitFac, this->absPrec);
    }
};

}
