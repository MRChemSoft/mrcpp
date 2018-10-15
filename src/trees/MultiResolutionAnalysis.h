#pragma once

#include "BoundingBox.h"
#include "core/ScalingBasis.h"
#include "core/MWFilter.h"

#include "mrcpp_declarations.h"

namespace mrcpp {

template<int D>
class MultiResolutionAnalysis final {
public:
    MultiResolutionAnalysis(const BoundingBox<D> &bb, const ScalingBasis &sb, int depth = MaxDepth);
    MultiResolutionAnalysis(const MultiResolutionAnalysis<D> &mra);
    MultiResolutionAnalysis &operator=(const MultiResolutionAnalysis &mra) = delete;

    int getOrder() const { return this->basis.getScalingOrder(); }
    int getMaxDepth() const { return this->maxDepth; }
    int getMaxScale() const { return this->world.getScale() + this->maxDepth; }

    const MWFilter &getFilter() const { return *this->filter; }
    const ScalingBasis &getScalingBasis() const { return this->basis; }
    const BoundingBox<D> &getWorldBox() const { return this->world; }

    MultiResolutionAnalysis<1> getKernelMRA() const;
    MultiResolutionAnalysis<2> getOperatorMRA() const;

    bool operator==(const MultiResolutionAnalysis<D> &mra) const;
    bool operator!=(const MultiResolutionAnalysis<D> &mra) const;

    void print() const;

protected:
    const int maxDepth;
    const ScalingBasis basis;
    const BoundingBox<D> world;
    MWFilter *filter;

    void setupFilter();
};

}
