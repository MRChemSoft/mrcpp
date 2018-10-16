#pragma once

#include "TreeAdaptor.h"
#include "constants.h"

/** Builds an OperatorTree with known band width (e.g. derivative and identity).
  * Assumes translational invariant and symmetric (in x - y) operator and keeps
  * only lower row of nodes (lx = 0). */

namespace mrcpp {

class BandWidthAdaptor final : public TreeAdaptor<2> {
public:
    BandWidthAdaptor(int bw, int ms) : TreeAdaptor<2>(ms), bandWidth(bw) { }

protected:
    const int bandWidth;

    bool splitNode(const MWNode<2> &node) const {
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

}
