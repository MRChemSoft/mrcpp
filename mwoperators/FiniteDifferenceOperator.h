#ifndef FINITEDIFFERENCEOPERATOR_H
#define FINITEDIFFERENCEOPERATOR_H

#include <Eigen/Core>
#include <string>

template<int D> class FunctionTree;

template<int D>
class FiniteDifferenceOperator {
public:
    FiniteDifferenceOperator(int m, int n);
    virtual ~FiniteDifferenceOperator() { }

    void operator()(FunctionTree<D> &out, FunctionTree<D> &inp);
protected:
    const int diff_order;
    const int approx_order;

    Eigen::VectorXd inp_data;
    Eigen::VectorXd out_data;

    Eigen::MatrixXd inp_pts;
    Eigen::MatrixXd out_pts;

    void compute();
    double computePoint(int M, int N, double x0);

    void plot(const std::string &file, Eigen::MatrixXd &pts, Eigen::VectorXd &vals);

    void setupQuadPoints(FunctionTree<D> &tree, Eigen::MatrixXd &pts);
    void compressTreeData(FunctionTree<D> &tree, Eigen::VectorXd &data);
    void expandTreeData(FunctionTree<D> &tree, Eigen::VectorXd &data);
};

#endif // FINITEDIFFERENCEOPERATOR_H
