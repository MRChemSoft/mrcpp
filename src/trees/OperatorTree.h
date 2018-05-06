#pragma once

#include "MWTree.h"

namespace mrcpp {

class OperatorTree: public MWTree<2> {
public:
    OperatorTree(const MultiResolutionAnalysis<2> &mra,
                 double np);
    virtual ~OperatorTree();

    double getNormPrecision() const { return this->normPrec; }

    void calcBandWidth(double prec = -1.0);
    void clearBandWidth();

    void setupOperNodeCache();
    void clearOperNodeCache();

    BandWidth &getBandWidth() { return *this->bandWidth; }
    const BandWidth &getBandWidth() const { return *this->bandWidth; }

    OperatorNode &getNode(int n, int l) { return *nodePtrAccess[n][l]; }
    const OperatorNode &getNode(int n, int l) const { return *nodePtrAccess[n][l]; }

    virtual void mwTransformDown(bool overwrite);
    virtual void mwTransformUp();

protected:
    const double normPrec;
    BandWidth *bandWidth;
    OperatorNode ***nodePtrStore;  ///< Avoids tree lookups
    OperatorNode ***nodePtrAccess; ///< Center (l=0) of node list

    void getMaxTranslations(Eigen::VectorXi &maxTransl);

    std::ostream& print(std::ostream &o);
};

}
