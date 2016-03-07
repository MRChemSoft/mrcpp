#ifndef MULTIRESOLUTIONANALYSIS_H
#define MULTIRESOLUTIONANALYSIS_H

#include "BoundingBox.h"
#include "ScalingBasis.h"
#include "FilterCache.h"

template<int D>
class MultiResolutionAnalysis {
public:
    MultiResolutionAnalysis(const MultiResolutionAnalysis<D> &mra)
            : world(mra.world),
              basis(mra.basis) {
        setupFilter();
    }
    MultiResolutionAnalysis(const BoundingBox<D> &bb, const ScalingBasis &sb)
            : world(bb),
              basis(sb) {
        setupFilter();
    }
    virtual ~MultiResolutionAnalysis() { }


    int getOrder() const { return this->basis.getScalingOrder(); }
    const MWFilter &getFilter() const { return *this->filter; }
    const ScalingBasis &getScalingBasis() const { return this->basis; }
    const BoundingBox<D> &getWorldBox() const { return this->world; }

    bool operator==(const MultiResolutionAnalysis<D> &mra) const {
        if (this->basis != mra.basis) return false;
        if (this->world != mra.world) return false;
        return true;
    }
    bool operator!=(const MultiResolutionAnalysis<D> &mra) const {
        if (this->basis != mra.basis) return true;
        if (this->world != mra.world) return true;
        return false;
    }

    /*
    template<int T>
    friend std::ostream& operator<<(std::ostream &o,
                                    const MultiResolutionAnalysis<T> &mra) {
        o << std::endl;
        o << "***************** MultiResolution Analysis *****************";
        o << std::endl;
        o << std::endl << mra.basis;
        o << std::endl;
        o << std::endl << mra.world;
        o << std::endl;
        o << "************************************************************";
        o << std::endl;
        return o;
    }
    */
protected:
    const ScalingBasis basis;
    const BoundingBox<D> world;
    MWFilter *filter;

    void setupFilter() {
        getLegendreFilterCache(lfilters);
        getInterpolatingFilterCache(ifilters);
        int k = this->basis.getScalingOrder();
        int type = this->basis.getScalingType();
        switch (type) {
        case Legendre:
            this->filter = &lfilters.get(k);
            break;
        case Interpol:
            this->filter = &ifilters.get(k);
            break;
        default:
            MSG_ERROR("Invalid scaling basis selected.")
        }
    }
};

#endif // MULTIRESOLUTIONANALYSIS_H
