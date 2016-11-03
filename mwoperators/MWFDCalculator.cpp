#include "MWFDCalculator.h"
#include "FunctionTree.h"
#include "FunctionNode.h"

using namespace std;
using namespace Eigen;

template<int D>
MWFDCalculator<D>::MWFDCalculator(int dir, int k, int n,
                                  FunctionTree<D> &f)
        : apply_dir(dir),
          diff_order(k),
          approx_order(n),
          f_tree(&f) {
}

template<int D>
void MWFDCalculator<D>::calcNode(MWNode<D> &node) {
    FunctionNode<D> &g_node = static_cast<FunctionNode<D> &>(node);
    const NodeIndex<D> &g_idx = g_node.getNodeIndex();

    MatrixXd input_data = setupInputData(g_idx);
    MatrixXd output_data = setupOutputData(g_node);
    compute(output_data, input_data);
    g_node.setValues(output_data.col(D));
}

template<int D>
MatrixXd MWFDCalculator<D>::setupInputData(const NodeIndex<D> &idx) {
    if (D != 1) NOT_IMPLEMENTED_ABORT;

    // Index for neighboring nodes in direction of application
    NodeIndex<D> idx_m1(idx);
    NodeIndex<D> idx_p1(idx);
    idx_m1.getTranslation()[this->apply_dir]--;
    idx_p1.getTranslation()[this->apply_dir]++;

    // Returns 0 if we go beyond world boundary
    int i_m1 = this->f_tree->getRootIndex(idx_m1);
    int i_0 = this->f_tree->getRootIndex(idx);
    int i_p1 = this->f_tree->getRootIndex(idx_p1);

    // Temp storage for quad pts and values
    MatrixXd pts;
    VectorXd vals;

    // Total number of quad pts (actually on children nodes)
    int n_coefs = this->f_tree->getTDim()*this->f_tree->getKp1_d();
    MatrixXd input_data = MatrixXd::Zero(3*n_coefs, D+1);

    if (i_0 < 0) MSG_FATAL("Node should exist!");
    MWNode<D> &mw_0 = this->f_tree->getNode(idx);
    FunctionNode<D> &f_0 = static_cast<FunctionNode<D> &>(mw_0);
    f_0.getPrimitiveChildPts(pts);
    f_0.getValues(vals);
    input_data.block(n_coefs,0,n_coefs,D) = pts;
    input_data.block(n_coefs,D,n_coefs,1) = vals;

    // Creating ghost points at boundaries (values will be zero)
    double n_size = pow(2.0, -idx.getScale());
    if (i_m1 < 0 or i_p1 < 0) {
        MatrixXd shift(n_coefs, D);
        shift.setConstant(n_size);
        if (i_m1 < 0) input_data.block(0*n_coefs,0,n_coefs,D) = pts - shift;
        if (i_p1 < 0) input_data.block(2*n_coefs,0,n_coefs,D) = pts + shift;
    }

    if (i_m1 >= 0) {
        // Fetch neighboring leaf node (only works for 1D)
        int fetch_scale = -1;
        double r_m1 = n_size*idx.getTranslation(0) - MachineZero;
        MWNode<D> &mw_m1 = this->f_tree->getNode(&r_m1, fetch_scale);
        FunctionNode<D> &f_m1 = static_cast<FunctionNode<D> &>(mw_m1);
        f_m1.getPrimitiveChildPts(pts);
        f_m1.getValues(vals);
        input_data.block(0,0,n_coefs,D) = pts;
        input_data.block(0,D,n_coefs,1) = vals;
    }
    if (i_p1 >= 0) {
        // Fetch neighboring leaf node (only works for 1D)
        int fetch_scale = -1;
        double r_p1 = n_size*idx_p1.getTranslation(0) + MachineZero;
        MWNode<D> &mw_p1 = this->f_tree->getNode(&r_p1, fetch_scale);
        FunctionNode<D> &f_p1 = static_cast<FunctionNode<D> &>(mw_p1);
        f_p1.getPrimitiveChildPts(pts);
        f_p1.getValues(vals);
        input_data.block(2*n_coefs,0,n_coefs,D) = pts;
        input_data.block(2*n_coefs,D,n_coefs,1) = vals;
    }
    return input_data;
}

template<int D>
MatrixXd MWFDCalculator<D>::setupOutputData(FunctionNode<D> &node) {
    int n_coefs = this->f_tree->getTDim()*this->f_tree->getKp1_d();
    MatrixXd output_data = MatrixXd::Zero(n_coefs, D+1);
    MatrixXd pts;
    node.getPrimitiveChildPts(pts);
    output_data.block(0,0,n_coefs,D) = pts;
    return output_data;
}

template<int D>
void MWFDCalculator<D>::compute(MatrixXd &output_data,
                                            MatrixXd &input_data) {
    int m = this->diff_order;
    int n = this->approx_order;
    int nOut = output_data.rows();
    for (int x = 0; x < nOut; x++) {
        double x0 = output_data(x, this->apply_dir);
        output_data(x, D) = computePoint(m, n, x0, input_data);
    }
}

template<int D>
double MWFDCalculator<D>::computePoint(int M,
                                                   int N,
                                                   double x0,
                                                   MatrixXd &input_data) {
    int d = this->apply_dir;
    int nInp = input_data.rows();
    int i0 = -1;
    for (int i = 0; i < nInp; i++) {
        if (fabs(x0 - input_data(i,d)) < MachineZero) i0 = i;
    }
    if (i0 < 0) MSG_FATAL("x0 did not coincide with an input point");

    int imin = N;
    int imax = nInp - 1 - N;

    if (i0 < imin) i0 = imin;
    if (i0 > imax) i0 = imax;

    VectorXd alpha = VectorXd::Zero(2*N+1);
    alpha(0) = input_data(i0, d);
    for (int nu = 1; nu <= N; nu++) {
        alpha(2*nu - 1) = input_data(i0 + nu, d);
        alpha(2*nu    ) = input_data(i0 - nu, d);
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
    double result = last_row(0)*input_data(i0, D);
    for (int nu = 1; nu <= N; nu++) {
        result += last_row(2*nu-1)*input_data(i0 + nu, D);
        result += last_row(2*nu  )*input_data(i0 - nu, D);
    }
    for (int m = 0; m <= M; m++) {
        delete delta[m];
    }
    delete[] delta;
    return result;
}

template class MWFDCalculator<1>;
template class MWFDCalculator<2>;
template class MWFDCalculator<3>;
