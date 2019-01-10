#pragma once

#include "TreeCalculator.h"
#include "operators/OperatorStatistics.h"

namespace mrcpp {

template<int D>
class DerivativeCalculator final : public TreeCalculator<D> {
public:
    DerivativeCalculator(int dir, DerivativeOperator<D> &o, FunctionTree<D> &f);
    ~DerivativeCalculator();

    MWNodeVector<D>* getInitialWorkVector(MWTree<D> &tree) const;

private:
    int applyDir;
    FunctionTree<D> *fTree;
    DerivativeOperator<D> *oper;

    std::vector<Timer> band_t;
    std::vector<Timer> calc_t;
    std::vector<Timer> norm_t;
    OperatorStatistics<D> operStat;

    MWNodeVector<D> makeOperBand(const MWNode<D> &gNode);

    void initTimers();
    void clearTimers();
    void printTimers() const;

    void calcNode(MWNode<D> &node);
    void postProcess() {
        printTimers();
        clearTimers();
        initTimers();
    }

    void applyOperator(OperatorState<D> &os);
    void tensorApplyOperComp(OperatorState<D> &os);
};

}
