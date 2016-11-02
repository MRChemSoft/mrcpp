#ifndef MWFDCALCULATOR_H
#define MWFDCALCULATOR_H

#include <Eigen/Core>

#include "TreeCalculator.h"

template<int D>
class MWFDCalculator : public TreeCalculator<D> {
public:
    MWFDCalculator(int dir, int k, int n, FunctionTree<D> &f);
protected:
    const int apply_dir;
    const int diff_order;
    const int approx_order;
    FunctionTree<D> *f_tree;

    virtual void calcNode(MWNode<D> &node);

    Eigen::MatrixXd setupInputData(const NodeIndex<D> &idx);
    Eigen::MatrixXd setupOutputData(FunctionNode<D> &node);
    void compute(Eigen::MatrixXd &output_data, Eigen::MatrixXd &input_data);
    double computePoint(int M, int N, double x0, Eigen::MatrixXd &input_data);
};

#endif // MWFDCALCULATOR_H
