#ifndef OPERATORADAPTOR_H
#define OPERATORADAPTOR_H

#include "WaveletAdaptor.h"

class OperatorAdaptor : public WaveletAdaptor<2> {
public:
    OperatorAdaptor(double p = -1.0, int n = MaxScale, bool a = false)
            : WaveletAdaptor<2>(p, n, a) { }
    virtual ~OperatorAdaptor() { }

protected:
    virtual bool splitNode(const MWNode<2> &node) const {
        int chkCompNorm = 0;
        for (int i = 1; i < 4; i++) {
            if (node.getComponentNorm(i) > 0.0) {
                chkCompNorm = 1;
            }
        }

        const int *l = node.getTranslation();
        int chkTransl = (l[0] == 0 or l[1] == 0);

        int scale = node.getScale();
        int chkScale = ((scale + 1) < this->maxScale);

        int split = 1;
        split *= chkTransl;
        split *= chkScale;
        split *= chkCompNorm;

        return split;
    }
};

#endif // OPERATORADAPTOR_H
