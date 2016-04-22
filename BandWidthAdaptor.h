#ifndef BANDWIDTHADAPTOR_H
#define BANDWIDTHADAPTOR_H

#include "TreeAdaptor.h"

/** Builds an OperatorTree with known band width (e.g. derivative and identity).
  * Assumes translational invariant and symmetric (in x - y) operator and keeps
  * only lower row of nodes (lx = 0). */

class BandWidthAdaptor : public TreeAdaptor<2> {
public:
    BandWidthAdaptor(int bw) : bandWidth(bw) { }
    BandWidthAdaptor(const BandWidthAdaptor &a) : bandWidth(a.bandWidth) { }
    virtual ~BandWidthAdaptor() { }
    virtual TreeAdaptor<2> *copy() const { return new BandWidthAdaptor(*this); }

protected:
    const int bandWidth;

    virtual bool splitNode(const MWNode<2> &node) const {
        const NodeIndex<2> &idx = node.getNodeIndex();
        int lx = idx.getTranslation(0);
        int ly = idx.getTranslation(1);
        int dl = abs(lx - ly);
        // Within band width on NEXT scale
        if ((lx == 0) and (2*dl <= this->bandWidth)) {
            return true;
        }
        return false;
    }
};

#endif // BANDWIDTHADAPTOR_H
