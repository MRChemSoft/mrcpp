#pragma once

#include "TreeCalculator.h"
#include "operators/OperatorStatistics.h"

#include "mrcpp_declarations.h"

namespace mrcpp {

template<int D>
class ConvolutionCalculator final : public TreeCalculator<D> {
public:
    ConvolutionCalculator(double p, ConvolutionOperator<D> &o, FunctionTree<D> &f, int depth = MaxDepth);
    ~ConvolutionCalculator();

    MWNodeVector<D>* getInitialWorkVector(MWTree<D> &tree) const;

private:
    int maxDepth;
    double prec;
    ConvolutionOperator<D> *oper;
    FunctionTree<D> *fTree;
    std::vector<Timer *> band_t;
    std::vector<Timer *> calc_t;
    std::vector<Timer *> norm_t;

    OperatorStatistics<D> operStat;
    std::vector<Eigen::MatrixXi *> bandSizes;

    static const int nComp = (1 << D);
    static const int nComp2 = (1 << D) * (1 << D);

    MWNodeVector<D>* makeOperBand(const MWNode<D> &gNode);
    void fillOperBand(MWNodeVector<D> *band, NodeIndex<D> &idx, const int *nbox, int dim);

    void initTimers();
    void clearTimers();
    void printTimers() const;

    void initBandSizes();
    int getBandSizeFactor(int i, int depth,const OperatorState<D> &os) const;
    void calcBandSizeFactor(Eigen::MatrixXi &bs, int depth, const BandWidth &bw);

    void calcNode(MWNode<D> &node);
    void postProcess() {
        printTimers();
        clearTimers();
        initTimers();
    }

    void applyOperComp(OperatorState<D> &os);
    void applyOperator(OperatorState<D> &os);
    void tensorApplyOperComp(OperatorState<D> &os);
};

}
