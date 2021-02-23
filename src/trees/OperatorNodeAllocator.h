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

/**
 *
 *
 *  \date Jul, 2016
 *  \author Peter Wind <peter.wind@uit.no> \n
 *  CTCC, University of Tromsø
 *
 */

#pragma once

#include "NodeAllocator.h"

namespace mrcpp {

class OperatorNodeAllocator final : public NodeAllocator<2> {
public:
    OperatorNodeAllocator(OperatorTree *tree);
    OperatorNodeAllocator(const OperatorNodeAllocator &tree) = delete;
    OperatorNodeAllocator &operator=(const OperatorNodeAllocator &tree) = delete;
    ~OperatorNodeAllocator() override;

    int alloc(int nAlloc, bool coeff = true);
    void dealloc(int serialIx) override;

    int getNChunks() const override;
    int getNodeChunkSize() const;
    int getCoeffChunkSize() const;

    double * getCoef_p(int sIdx);
    OperatorNode * getNode_p(int sIdx);

protected:
    OperatorNode *lastNode{nullptr};     // pointer to the last active node
    std::vector<OperatorNode *> nodeChunks{};

    double * getCoefNoLock(int sIdx);
    OperatorNode * getNodeNoLock(int sIdx);

    void appendChunk(bool coeff = true);
};

} // namespace mrcpp
