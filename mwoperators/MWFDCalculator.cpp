#include "MWFDCalculator.h"
#include "FunctionTree.h"
#include "FunctionNode.h"
#include "MathUtils.h"

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

    VectorXd inp_pts;
    VectorXd out_pts;
    MatrixXd inp_vals;
    MatrixXd out_vals;

    fetchInputData(g_idx, inp_pts, inp_vals);
    fetchOutputPoints(g_node, out_pts);
    computeBlock(out_pts, out_vals, inp_pts, inp_vals);
    pushOutputValues(g_node, out_vals);
}

template<int D>
void MWFDCalculator<D>::fetchInputData(const NodeIndex<D> &idx,
                                       VectorXd &pts,
                                       MatrixXd &vals) {
    // Index for neighboring nodes in direction of application
    NodeIndex<D> idx_m1(idx);
    NodeIndex<D> idx_p1(idx);
    idx_m1.getTranslation()[this->apply_dir]--;
    idx_p1.getTranslation()[this->apply_dir]++;

    // Returns 0 if we go beyond world boundary
    int i_m1 = this->f_tree->getRootIndex(idx_m1);
    int i_0 = this->f_tree->getRootIndex(idx);
    int i_p1 = this->f_tree->getRootIndex(idx_p1);

    // Total number of quad pts (actually on children nodes, hence 2*)
    int tkp1 = 2*this->f_tree->getKp1();
    int tkp1_dm1 = MathUtils::ipow(tkp1, D-1);
    pts = VectorXd::Zero(3*tkp1);
    vals = MatrixXd::Zero(tkp1_dm1, 3*tkp1);

    if (i_0 >= 0) {
        VectorXd node_pts;
        MatrixXd node_vals = MatrixXd::Zero(tkp1_dm1, tkp1);
        MWNode<D> &mw_0 = this->f_tree->getNode(idx);
        FunctionNode<D> &f_0 = static_cast<FunctionNode<D> &>(mw_0);
        getNodeData(f_0, node_pts, node_vals);
        pts.segment(1*tkp1, tkp1) = node_pts;
        vals.block(0, 1*tkp1, tkp1_dm1, tkp1) = node_vals;
    } else {
        MSG_FATAL("Node should exist!");
    }

    // Creating ghost points at boundaries (values will be zero)
    double n_size = pow(2.0, -idx.getScale());
    VectorXd shift(tkp1);
    shift.setConstant(n_size);

    if (i_m1 >= 0) {
        VectorXd node_pts;
        MatrixXd node_vals = MatrixXd::Zero(tkp1_dm1, tkp1);
        MWNode<D> &mw_m1 = this->f_tree->getNode(idx_m1);
        FunctionNode<D> &f_m1 = static_cast<FunctionNode<D> &>(mw_m1);
        getNodeData(f_m1, node_pts, node_vals);
        pts.segment(0*tkp1, tkp1) = node_pts;
        vals.block(0, 0*tkp1, tkp1_dm1, tkp1) = node_vals;
    } else {
        pts.segment(0, tkp1) = pts.segment(tkp1, tkp1) - shift;
    }

    if (i_p1 >= 0) {
        VectorXd node_pts;
        MatrixXd node_vals = MatrixXd::Zero(tkp1_dm1, tkp1);
        MWNode<D> &mw_p1 = this->f_tree->getNode(idx_p1);
        FunctionNode<D> &f_p1 = static_cast<FunctionNode<D> &>(mw_p1);
        getNodeData(f_p1, node_pts, node_vals);
        pts.segment(2*tkp1, tkp1) = node_pts;
        vals.block(0, 2*tkp1, tkp1_dm1, tkp1) = node_vals;
    } else {
        pts.segment(2*tkp1,tkp1) = pts.segment(tkp1, tkp1) + shift;
    }
}

template<int D>
void MWFDCalculator<D>::getNodeData(FunctionNode<D> &node,
                                    VectorXd &pts,
                                    MatrixXd &vals) {
    NOT_IMPLEMENTED_ABORT;
}

template<>
void MWFDCalculator<3>::getNodeData(FunctionNode<3> &node,
                                    VectorXd &pts,
                                    MatrixXd &vals) {
    int dir = this->apply_dir;
    int kp1 = node.getKp1();
    int kp1_dm1 = MathUtils::ipow(kp1, 2);
    int kp1_d = node.getKp1_d();
    int tDim = node.getTDim();

    MatrixXd prim_pts;
    node.getPrimitiveChildPts(prim_pts);
    pts = prim_pts.row(dir);

    VectorXd unsort_vals;
    node.getValues(unsort_vals);

    int X,Y,Z;
    MatrixXd **sort_vals = new MatrixXd*[tDim];
    for (int t = 0; t < tDim; t++) {
        sort_vals[t] = new MatrixXd;
        const VectorXd &inp = unsort_vals.segment(t*kp1_d, kp1_d);
        MatrixXd &out = *sort_vals[t];

        out = MatrixXd::Zero(kp1_dm1, kp1);
        int n = 0;
        for (int x = 0; x < kp1; x++) {
            for (int y = 0; y < kp1; y++) {
                for (int z = 0; z < kp1; z++) {
                    if (dir == 0) { X = x; Y = y; Z = z; }
                    if (dir == 1) { X = y; Y = z; Z = x; }
                    if (dir == 2) { X = z; Y = x; Z = y; }
                    out(n,z) = inp(X*kp1*kp1 + Y*kp1 + Z);
                }
                n++;
            }
        }
    }
    for (int z = 0; z < 2; z++) {
        for (int y = 0; y < 2; y++) {
            for (int x = 0; x < 2; x++) {
                if (dir == 0) { X = x; Y = y; Z = z; }
                if (dir == 1) { X = y; Y = z; Z = x; }
                if (dir == 2) { X = z; Y = x; Z = y; }
                vals.block((2*z+y)*kp1_dm1, x*kp1, kp1_dm1, kp1) = *sort_vals[4*Z+2*Y+X];
            }
        }
    }

    for (int t = 0; t < tDim; t++) {
        delete sort_vals[t];
    }
    delete[] sort_vals;
}

template<int D>
void MWFDCalculator<D>::fetchOutputPoints(const MWNode<D> &node,
                                          VectorXd &pts) {
    MatrixXd prim_pts;
    node.getPrimitiveChildPts(prim_pts);
    pts = prim_pts.row(this->apply_dir);
}

template<int D>
void MWFDCalculator<D>::pushOutputValues(FunctionNode<D> &node,
                                         MatrixXd &vals) {
    NOT_IMPLEMENTED_ABORT;
}

template<>
void MWFDCalculator<3>::pushOutputValues(FunctionNode<3> &node,
                                         MatrixXd &vals) {
    int dir = this->apply_dir;
    int kp1 = node.getKp1();
    int kp1_dm1 = MathUtils::ipow(kp1, 2);
    int kp1_d = node.getKp1_d();
    int tDim = node.getTDim();

    MatrixXd **sort_vals = new MatrixXd*[tDim];
    for (int t = 0; t < tDim; t++) {
        sort_vals[t] = new MatrixXd;
    }

    int X,Y,Z;
    for (int z = 0; z < 2; z++) {
        for (int y = 0; y < 2; y++) {
            for (int x = 0; x < 2; x++) {
                if (dir == 0) { X = x; Y = y; Z = z; }
                if (dir == 1) { Y = x; Z = y; X = z; }
                if (dir == 2) { Z = x; X = y; Y = z; }
                *sort_vals[4*Z+2*Y+X] = vals.block((2*z+y)*kp1_dm1, x*kp1, kp1_dm1, kp1);
            }
        }
    }

    VectorXd unsort_vals = VectorXd::Zero(tDim*kp1_d);
    for (int t = 0; t < tDim; t++) {
        MatrixXd &inp = *sort_vals[t];
        VectorXd out = VectorXd::Zero(kp1_d);
        int n = 0;
        for (int x = 0; x < kp1; x++) {
            for (int y = 0; y < kp1; y++) {
                for (int z = 0; z < kp1; z++) {
                    if (dir == 0) { X = x; Y = y; Z = z; }
                    if (dir == 1) { X = y; Y = z; Z = x; }
                    if (dir == 2) { X = z; Y = x; Z = y; }
                    out(X*kp1*kp1 + Y*kp1 + Z) = inp(n,z);
                }
                n++;
            }
        }
        unsort_vals.segment(t*kp1_d, kp1_d) = out;
    }
    for (int t = 0; t < tDim; t++) {
        delete sort_vals[t];
    }
    delete[] sort_vals;

    node.setCoefBlock(0, tDim*kp1_d, unsort_vals.data());
    node.cvTransform(Backward);
    node.mwTransform(Compression);
    node.setHasCoefs();
    node.calcNorms();
}

template<int D>
void MWFDCalculator<D>::computeBlock(VectorXd &out_pts,
                                     MatrixXd &out_vals,
                                     VectorXd &inp_pts,
                                     MatrixXd &inp_vals) {
    out_vals = MatrixXd::Zero(inp_vals.rows(), out_pts.size());
    for (int i = 0; i < inp_vals.rows(); i++) {
        MatrixXd out_data = MatrixXd::Zero(out_pts.size(), 2);
        MatrixXd inp_data = MatrixXd::Zero(inp_pts.size(), 2);
        out_data.col(0) = out_pts;
        inp_data.col(0) = inp_pts;
        inp_data.col(1) = inp_vals.row(i);
        computeLine(out_data, inp_data);
        out_vals.row(i) = out_data.col(1);
    }
}

template<int D>
void MWFDCalculator<D>::computeLine(MatrixXd &output_data,
                                    MatrixXd &input_data) {
    int m = this->diff_order;
    int n = this->approx_order;
    int nOut = output_data.rows();
    for (int x = 0; x < nOut; x++) {
        double x0 = output_data(x, 0);
        output_data(x, 1) = computePoint(m, n, x0, input_data);
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
        if (fabs(x0 - input_data(i,0)) < MachineZero) i0 = i;
    }
    if (i0 < 0) MSG_FATAL("x0 did not coincide with an input point");

    int imin = N;
    int imax = nInp - 1 - N;

    if (i0 < imin) i0 = imin;
    if (i0 > imax) i0 = imax;

    VectorXd alpha = VectorXd::Zero(2*N+1);
    alpha(0) = input_data(i0, 0);
    for (int nu = 1; nu <= N; nu++) {
        alpha(2*nu - 1) = input_data(i0 + nu, 0);
        alpha(2*nu    ) = input_data(i0 - nu, 0);
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
    double result = last_row(0)*input_data(i0, 1);
    for (int nu = 1; nu <= N; nu++) {
        result += last_row(2*nu-1)*input_data(i0 + nu, 1);
        result += last_row(2*nu  )*input_data(i0 - nu, 1);
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
