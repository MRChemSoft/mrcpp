#ifndef OPERAPPLICATIONCALCULATOR_H
#define OPERAPPLICATIONCALCULATOR_H

#include <map>

#include "TreeCalculator.h"
#include "mwrepr_declarations.h"

template<int D>
class OperApplicationCalculator : public TreeCalculator<D> {
public:
    OperApplicationCalculator(OperatorTreeVector &o, FunctionTree<D> &f);
    virtual ~OperApplicationCalculator();

    virtual MWNodeVector* getInitialWorkVector(MWTree<D> &tree) const;

protected:
    OperatorTreeVector *oper;
    FunctionTree<D> *fTree;

    // Operator application statistics
    int nThreads;
    Eigen::Matrix<int, 8, 8> **compCount;
    Eigen::Vector2i **gNodeCount;
    Eigen::Vector2i **fNodeCount;
    Eigen::Vector2i totFCount;
    Eigen::Vector2i totGCount;
    int *genCount;
    int totAppNodes;
    int totGenAppNodes;
    std::vector<Eigen::MatrixXi *> bandSizes;

    static const int nComp = (1 << D);
    static const int nComp2 = (1 << D) * (1 << D);

    void initBandSizes();

    void initNodeCounters();
    void resetNodeCounters();
    void flushNodeCounters();
    void incrementNodeCounters(int ft, int gt, bool isGenNode);
    void incrementFNodeCounters(MWNode<D> &gNode);
    void incrementGNodeCounters(MWNode<D> &gNode);

    MWNodeVector* makeOperBand(const MWNode<D> &gNode);
    void fillOperBand(MWNodeVector *band, NodeIndex<D> &idx, const int *nbox, int dim);
    void freeOperBand(MWNodeVector *band, MWNode<D> &gNode);

    int getBandSizeFactor(int i, int depth,const OperatorState<D> &os) const;
    void calcBandSizeFactor(Eigen::MatrixXi &bs, int depth, const BandWidth &bw);

    virtual void calcNode(MWNode<D> &node);

    void applyOperComp(OperatorState<D> &os);
    void applyOperator(OperatorState<D> &os);
    void tensorApplyOperComp(OperatorState<D> &os);
};


#endif // OPERAPPLICATIONCALCULATOR_H
