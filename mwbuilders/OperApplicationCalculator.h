#ifndef OPERAPPLICATIONCALCULATOR_H
#define OPERAPPLICATIONCALCULATOR_H

#include "TreeCalculator.h"
#include "OperatorStatistics.h"
#include "mrcpp_declarations.h"

template<int D>
class OperApplicationCalculator : public TreeCalculator<D> {
public:
    OperApplicationCalculator(int dir,
                              double p,
                              OperatorTreeVector &o,
                              FunctionTree<D> &f,
                              int depth = MaxDepth);
    virtual ~OperApplicationCalculator();

    virtual MWNodeVector* getInitialWorkVector(MWTree<D> &tree) const;

protected:
    int applyDir;
    int maxDepth;
    double prec;
    OperatorTreeVector *oper;
    FunctionTree<D> *fTree;

    OperatorStatistics<D> operStat;
    std::vector<Eigen::MatrixXi *> bandSizes;

    static const int nComp = (1 << D);
    static const int nComp2 = (1 << D) * (1 << D);

    MWNodeVector* makeOperBand(const MWNode<D> &gNode);
    void fillOperBand(MWNodeVector *band, NodeIndex<D> &idx, const int *nbox, int dim);
    void freeOperBand(MWNodeVector *band, MWNode<D> &gNode);

    void initBandSizes();
    int getBandSizeFactor(int i, int depth,const OperatorState<D> &os) const;
    void calcBandSizeFactor(Eigen::MatrixXi &bs, int depth, const BandWidth &bw);

    virtual void calcNode(MWNode<D> &node);

    void applyOperComp(OperatorState<D> &os);
    void applyOperator(OperatorState<D> &os);
    void tensorApplyOperComp(OperatorState<D> &os);
};


#endif // OPERAPPLICATIONCALCULATOR_H
