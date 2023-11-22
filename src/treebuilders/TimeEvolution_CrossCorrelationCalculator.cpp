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
 * 
 */
void TimeEvolution_CrossCorrelationCalculator::calcNode(MWNode<2> &node)
{
    node.zeroCoefs();
    int type = node.getMWTree().getMRA().getScalingBasis().getScalingType();
    switch (type) {
        case Interpol: {
            MSG_ERROR("Invalid scaling type");
            break;
        }
        case Legendre: {
            applyCcc(node);
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
 *  @details This will ... (work in progress)
 * 
 * 
 * 
 * 
 */
//template <int T>
void TimeEvolution_CrossCorrelationCalculator::applyCcc(MWNode<2> &node)
{
    // The scale of J power integrals:
    //int scale = node.getScale() + 1;  //scale = n = (n - 1) + 1
    
    int t_dim = node.getTDim();       //t_dim = 4
    int kp1_d = node.getKp1_d();      //kp1_d = (k + 1)^2

    VectorXd vec_o = VectorXd::Zero(t_dim * kp1_d);
    const NodeIndex<2> &idx = node.getNodeIndex();
    
    for (int i = 0; i < t_dim; i++)
    {
        NodeIndex<2> l = idx.child(i);
        int l_b = l[1] - l[0];

        int vec_o_segment_index = 0;
        for( int p = 0; p <= node.getOrder(); p++ )
            for( int j = 0; j <= node.getOrder(); j++ )
            {
                //std::min(M, N)  could be used for breaking the following loop
                //this->cross_correlation->Matrix.size() should be big enough a priori
                for( int k = 0; 2*k + p + j < (*this->J_power_inetgarls)[l_b].size(); k++ )
                {
                    double J;
                    if( this->imaginary ) J = (*this->J_power_inetgarls)[l_b][2*k + p + j].imag();
                    else J = (*this->J_power_inetgarls)[l_b][2*k + p + j].real();
                    vec_o.segment(i * kp1_d, kp1_d)(vec_o_segment_index)
                    +=
                    J * cross_correlation->Matrix[k](p, j); //by default eigen library reads a transpose matrix from a file
                }
                vec_o_segment_index++;
            }
    }

    double *coefs = node.getCoefs();
    for (int i = 0; i < t_dim * kp1_d; i++) {
        //auto scaling_factor = node.getMWTree().getMRA().getWorldBox().getScalingFactor(0);
        coefs[i] = vec_o(i);
        //std::cout<< "coefs[i] = " << coefs[i] << std::endl;
    }
}

} // namespace mrcpp
