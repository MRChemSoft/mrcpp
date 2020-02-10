/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2020 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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

#pragma once

#include "TreeAdaptor.h"

namespace mrcpp {

template <int D> class WaveletAdaptor : public TreeAdaptor<D> {
public:
    WaveletAdaptor(double pr, int ms, bool ap = false, double sf = 1.0)
            : TreeAdaptor<D>(ms)
            , absPrec(ap)
            , prec(pr)
            , splitFac(sf) {}
    ~WaveletAdaptor() override = default;

    void setPrecTree(std::vector<FunctionTree<D> *> treevec) { this->precTrees = treevec; }
    void setMultiplicationSplit(bool multSplit) { this->multiplicationSplit = multSplit; };

protected:
    bool absPrec;
    double prec;
    double splitFac;
    std::vector<FunctionTree<D> *> precTrees;
    bool multiplicationSplit = false;

    bool splitNode(const MWNode<D> &node) const override {
      if(multiplicationSplit and this->precTrees.size()>0){
	auto multPrec = 1.0;
        if(this->precTrees.size()!=2)std::cout<<this->precTrees.size()<<" ERROR "<<std::endl;
        auto &pNode0 = precTrees[0]->getNode(node.getNodeIndex());
        auto &pNode1 = precTrees[1]->getNode(node.getNodeIndex());
        double maxW0=std::sqrt(pNode0.getMaxWSquareNorm());
        double maxW1=std::sqrt(pNode1.getMaxWSquareNorm());
        double maxS0=std::sqrt(pNode0.getMaxSquareNorm());
        double maxS1=std::sqrt(pNode1.getMaxSquareNorm());
        if(pNode0.isGenNode()){
            maxW0 = 0.0;
            maxS0 = std::sqrt(std::pow(2.0, D * node.getScale()) *pNode0.getSquareNorm());
        }
        if(pNode1.isGenNode()){
            maxW1 = 0.0;
            maxS1 = std::sqrt(std::pow(2.0, D * node.getScale()) *pNode1.getSquareNorm());
        }
        //The wavelet contribution (in the product of node0 and node1) can be approximated as
        multPrec =  maxW0*maxS1 +  maxW1*maxS0 + maxW0*maxW1 ;

        // Note: this never refine deeper than one scale more than input tree grids, because when wavelets are zero for both input trees, multPrec=0
        // In addition, we force not to refine deeper than input tree grids
	if (multPrec > this->prec and not (pNode0.isLeafNode() and pNode1.isLeafNode())) {
	  return true;
	} else {
	  return false;
	}
      }else{
        auto precNorm = 1.0;
        if(this->precTrees.size() >0) precNorm = 0.0;//initialize
        for (int i = 0; i < this->precTrees.size(); i++) {
            auto &pNode = precTrees[i]->getNode(node.getNodeIndex());
            auto n = node.getScale();
            if(not pNode.isGenNode() and pNode.getMaxSquareNorm()>0.0){
                precNorm = std::max(precNorm, std::sqrt(pNode.getMaxSquareNorm()));
            }else{
                precNorm = std::max(precNorm, std::sqrt(std::pow(2.0, D * n) * pNode.getSquareNorm()));
            }
	}
        return node.splitCheck(this->prec / precNorm, this->splitFac, this->absPrec);
      }
    }
};

} // namespace mrcpp
