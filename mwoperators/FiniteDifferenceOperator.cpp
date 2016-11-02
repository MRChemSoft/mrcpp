#include <fstream>

#include "FiniteDifferenceOperator.h"
#include "FunctionTree.h"
#include "FunctionNode.h"

using namespace std;
using namespace Eigen;

template<int D>
FiniteDifferenceOperator<D>::FiniteDifferenceOperator(int m, int n)
        : diff_order(m),
          approx_order(n) {
    if (this->approx_order < this->diff_order) MSG_FATAL("Invalid FD oper");
}

template<int D>
void FiniteDifferenceOperator<D>::operator()(FunctionTree<D> &out,
                                             FunctionTree<D> &inp) {
    if (D != 1) NOT_IMPLEMENTED_ABORT;
    setupQuadPoints(inp, this->inp_pts);
    setupQuadPoints(out, this->out_pts);
    compressTreeData(inp, this->inp_data);
    compute();
    expandTreeData(out, this->out_data);
}

template<int D>
void FiniteDifferenceOperator<D>::setupQuadPoints(FunctionTree<D> &tree,
                                                  MatrixXd &pts) {
    int nNodes = tree.getNEndNodes();
    int nCoefs = (1 << D)*tree.getKp1_d();
    pts = MatrixXd::Zero(nNodes*nCoefs, D);
    for (int n = 0; n < nNodes; n++) {
        FunctionNode<D> &node = tree.getEndFuncNode(n);
        MatrixXd x_i;
        node.getQuadraturePoints(x_i);
        pts.block(n*nCoefs, 0, nCoefs, D) = x_i;
    }
}

template<int D>
void FiniteDifferenceOperator<D>::compressTreeData(FunctionTree<D> &tree,
                                                   VectorXd &data) {
    int nNodes = tree.getNEndNodes();
    int nCoefs = (1 << D)*tree.getKp1_d();
    data = VectorXd::Zero(nNodes*nCoefs);
    for (int n = 0; n < nNodes; n++) {
        FunctionNode<D> &node = tree.getEndFuncNode(n);
        VectorXd c_i;
        node.getValues(c_i);
        data.segment(n*nCoefs, nCoefs) = c_i;
    }
}

template<int D>
void FiniteDifferenceOperator<D>::expandTreeData(FunctionTree<D> &tree,
                                                 VectorXd &data) {
    int nNodes = tree.getNEndNodes();
    int nCoefs = (1 << D)*tree.getKp1_d();
    if (data.size() != nNodes*nCoefs) MSG_FATAL("Non matching grids");
    for (int n = 0; n < nNodes; n++) {
        FunctionNode<D> &node = tree.getEndFuncNode(n);
        VectorXd c_i = data.segment(n*nCoefs, nCoefs);
        node.setValues(c_i);
    }
    tree.mwTransform(BottomUp);
    tree.calcSquareNorm();
}

template<int D>
double FiniteDifferenceOperator<D>::computePoint(int M, int N, double x0) {
    int nInp = this->inp_pts.size();
    int i0 = -1;
    for (int i = 0; i < nInp; i++) {
        if (fabs(x0 - this->inp_pts(i)) < MachineZero) i0 = i;
    }
    if (i0 < 0) MSG_FATAL("x0 did not coincide with an input point");

    int imin = N;
    int imax = this->inp_pts.size() - 1 - N;

    if (i0 < imin) i0 = imin;
    if (i0 > imax) i0 = imax;

    VectorXd alpha = VectorXd::Zero(2*N+1);
    alpha(0) = this->inp_pts(i0);
    for (int nu = 1; nu <= N; nu++) {
        alpha(2*nu - 1) = this->inp_pts(i0 + nu);
        alpha(2*nu    ) = this->inp_pts(i0 - nu);
    }

    MatrixXd **delta = new MatrixXd*[M+1];
    for (int m = 0; m <= M; m++) {
        delta[m] = new MatrixXd;
        *delta[m] = MatrixXd::Zero(2*N+1, 2*N+1);
    }

    (*delta[0])(0,0) = 1.0;
    double c1 = 1.0;
    for (int n = 1; n <= 2*N; n++) {
        double c2 = 1.0;
        for (int nu = 0; nu < n; nu++) {
            double c3 = alpha(n) - alpha(nu);
            c2 *= c3;
            if (n <= M) (*delta[n])(n-1,nu) = 0.0;
            (*delta[0])(n,nu) = (alpha(n) - x0)*(*delta[0])(n-1,nu)/c3;
            for (int m = 1; m <= min(n, M); m++) {
                (*delta[m])(n,nu) = ((alpha(n) - x0)*(*delta[m])(n-1,nu)
                                     -m*(*delta[m-1])(n-1,nu))/c3;
            }
        }
        (*delta[0])(n,n) = -(c1/c2)*(alpha(n-1) - x0)*(*delta[0])(n-1,n-1);
        for (int m = 1; m <= min(n, M); m++) {
            (*delta[m])(n,n) = (c1/c2)*(m*(*delta[m-1])(n-1,n-1)
                                      - (alpha(n-1) - x0)*(*delta[m])(n-1,n-1));
        }
        c1 = c2;
    }
    const MatrixXd &delta_M = (*delta[M]);
    const VectorXd &last_row = delta_M.row(delta_M.rows() - 1);
    double result = last_row(0)*this->inp_data(i0);
    for (int nu = 1; nu <= N; nu++) {
        result += last_row(2*nu-1)*this->inp_data(i0 + nu);
        result += last_row(2*nu  )*this->inp_data(i0 - nu);
    }
    for (int m = 0; m <= M; m++) {
        delete delta[m];
    }
    delete[] delta;
    return result;
}

template<int D>
void FiniteDifferenceOperator<D>::compute() {
    int m = this->diff_order;
    int n = this->approx_order;
    int nOut = this->out_pts.size();
    this->out_data = VectorXd::Zero(nOut);
    for (int x = 0; x < nOut; x++) {
        double x0 = this->out_pts(x);
        this->out_data(x) = computePoint(m, n, x0);
    }
}

template<int D>
void FiniteDifferenceOperator<D>::plot(const string &file,
                                       MatrixXd &pts,
                                       VectorXd &vals) {
    if (pts.size() != vals.size()) MSG_FATAL("Invalid arguments");

    ofstream of;
    of.open(file.c_str());
    if (of.bad()) MSG_FATAL("File error");

    for (int i = 0; i < pts.size(); i++) {
        of << scientific << setprecision(14);
        for (int d = 0; d < D; d++) {
            of << setw(25) << pts(i, d);
        }
        of << setw(25) << vals(i) << endl;
    }
    of.close();
}


template class FiniteDifferenceOperator<1>;
template class FiniteDifferenceOperator<2>;
template class FiniteDifferenceOperator<3>;
