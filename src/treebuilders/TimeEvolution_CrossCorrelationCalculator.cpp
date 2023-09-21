/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2021 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
 *
 * This file is part of MRCPP.
 *
 * MRCPP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRCPP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRCPP.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRCPP, see:
 * <https://mrcpp.readthedocs.io/>
 */

#include "TimeEvolution_CrossCorrelationCalculator.h"
#include "trees/FunctionTree.h"
#include "trees/MWNode.h"
#include "utils/Printer.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace mrcpp {

/** @param[in] node: ...
 *  @details This will ... (work in progress)
 * 
 * 
 * 
 */
void TimeEvolution_CrossCorrelationCalculator::calcNode(MWNode<2> &node)
{
    node.zeroCoefs();
    int type = node.getMWTree().getMRA().getScalingBasis().getScalingType();
    switch (type) {
        case Interpol: {
            getCrossCorrelationCache(Interpol, ccc);
            applyCcc(node, ccc);
            break;
        }
        case Legendre: {
            getCrossCorrelationCache(Legendre, ccc);
            applyCcc(node, ccc);
            break;
        }
        default:
            MSG_ERROR("Invalid scaling type");
            break;
    }
    node.mwTransform(Compression);
    node.setHasCoefs();
    node.calcNorms();
}


/** @param[in] node: ...
 *  @param[in] ccc: ...
 *  @details This will ... (work in progress)
 * 
 * 
 * 
 */
template <int T>
void TimeEvolution_CrossCorrelationCalculator::applyCcc(MWNode<2> &node, CrossCorrelationCache<T> &ccc)
{
    const MatrixXd &lMat = ccc.getLMatrix(node.getOrder());
    const MatrixXd &rMat = ccc.getRMatrix(node.getOrder());

    //std::cout << "Number of rows of lMat: " << lMat.rows() << std::endl;
    //std::cout << "Number of columns of lMat: " << lMat.cols() << std::endl;
    
    //std::cout << "lMat = " << std::endl << lMat << std::endl;
    //std::cout << "rMat = " << std::endl << std::setprecision(2) << rMat << std::endl;

    //std::cout << "The last column of rMat = " << std::endl;
    //for (int i = 0; i < rMat.rows(); i++)
    //{
    //    std::cout << rMat(i, rMat.cols() - 1) << std::endl;
    //}

    int scale = node.getScale() + 1;
    int t_dim = node.getTDim();
    int kp1_d = node.getKp1_d();

    std::cout<< "scale = n = (n - 1) + 1 = " << scale << std::endl;
    std::cout<< "t_dim = " << t_dim << std::endl;
    std::cout<< "kp1_d = " << kp1_d << std::endl;
    
    VectorXd vec_o = VectorXd::Zero(t_dim * kp1_d);
    const NodeIndex<2> &idx = node.getNodeIndex();
    
    std::cout<< "idx = " << idx << std::endl;
    //std::cout<< "t_dim * kp1_d = " << t_dim * kp1_d << std::endl;
    //std::cout<< "vec_o = " << vec_o << std::endl;

    for (int i = 0; i < t_dim; i++)
    {
        NodeIndex<2> l = idx.child(i);
        std::cout<< "idx.child(i) = " << idx.child(i) << std::endl;
        int l_a = l[1] - l[0] - 1;
        int l_b = l[1] - l[0];
        std::cout<< "l[0] = " << l[0] << std::endl;
        std::cout<< "l[1] = " << l[1] << std::endl;
        std::cout<< "l_a = " << l_a << std::endl;
        std::cout<< "l_b = " << l_b << std::endl;

        NodeIndex<1> idx_a(scale, {l_a});
        NodeIndex<1> idx_b(scale, {l_b});
        //std::cout<< "idx_a = " << idx_a << std::endl;
        //std::cout<< "idx_b = " << idx_b << std::endl;

        const MWNode<1> &node_a = this->kernel->getNode(idx_a);
        const MWNode<1> &node_b = this->kernel->getNode(idx_b);

        VectorXd vec_a;
        VectorXd vec_b;
        node_a.getCoefs(vec_a);
        node_b.getCoefs(vec_b);
        
        //std::cout<< "vec_a = " << vec_a << std::endl;
        //std::cout<< "vec_b = " << vec_b << std::endl;

        const VectorXd &seg_a = vec_a.segment(0, node_a.getKp1_d());
        const VectorXd &seg_b = vec_b.segment(0, node_b.getKp1_d());
        
        //std::cout<< "seg_a = " << seg_a << std::endl;
        //std::cout << "seg_b = " << std::endl << seg_b << std::endl;
        //std::cout << "node_b.getKp1_d() = " << node_b.getKp1_d() << std::endl;

        vec_o.segment(i * kp1_d, kp1_d) = (lMat * seg_a + rMat * seg_b);

        //std::cout << "seg_b.size() = " << seg_b.size() << std::endl;
        //std::cout << "Number of columns of rMat: " << rMat.cols() << std::endl;
        //std::cout<< "vec_o.segment(i * kp1_d, kp1_d) = " << vec_o.segment(i * kp1_d, kp1_d) << std::endl;
        
        //std::cout<< "vec_o = " << vec_o << std::endl;
        
    }
    double *coefs = node.getCoefs();
    double two_n = std::pow(2.0, -scale / 2.0);
    for (int i = 0; i < t_dim * kp1_d; i++) {
        auto scaling_factor = node.getMWTree().getMRA().getWorldBox().getScalingFactor(0);
        // This is only implemented for unifrom scaling factors
        // hence the zero TODO: make it work for non-unifrom scaling
        coefs[i] = std::sqrt(scaling_factor) * two_n * vec_o(i);

        //std::cout<< "coefs[i] = " << coefs[i] << std::endl;
    }

    //Can we create a complex operator?
    //std::complex<double> complexNum0(3.0, 4.0);  // Creates a complex number 3 + 4i
    //coefs[0] = complexNum0.real();

}

} // namespace mrcpp
