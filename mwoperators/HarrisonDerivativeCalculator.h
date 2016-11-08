#ifndef HARRISONDERIVATIVECALCULATOR_H
#define HARRISONDERIVATIVECALCULATOR_H

#include <Eigen/Core>

#include "TreeCalculator.h"

template<int D>
class HarrisonDerivativeCalculator : public TreeCalculator<D> {
public:
    HarrisonDerivativeCalculator(int dir, FunctionTree<D> &f);
    virtual ~HarrisonDerivativeCalculator() { }

protected:
    const int apply_dir;
    FunctionTree<D> *f_tree;

    MWNode<D> *f_m1;
    MWNode<D> *f_0;
    MWNode<D> *f_p1;

    Eigen::MatrixXd oper_m1;
    Eigen::MatrixXd oper_0;
    Eigen::MatrixXd oper_p1;

    virtual void calcNode(MWNode<D> &node);

    void initialize(int kp1);
    void fetchInputNodes(const NodeIndex<D> &idx_0);
    void clearInputNodes();
};

#endif // HARRISONDERIVATIVECALCULATOR_H
