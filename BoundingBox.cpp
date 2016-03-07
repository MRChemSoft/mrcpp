/**
 *
 *
 *          CTCC, University of Tromsø
 *
 */

#include "MathUtils.h"
#include "constants.h"
#include "BoundingBox.h"

using namespace std;

template<int D>
BoundingBox<D>::BoundingBox() {
    setNBoxes(0);
    setDerivedParameters();
}

template<int D>
BoundingBox<D>::BoundingBox(const NodeIndex<D> &idx, const int *nb)
            : cornerIndex(idx) {
    setNBoxes(nb);
    setDerivedParameters();
}

template<int D>
BoundingBox<D>::BoundingBox(const BoundingBox<D> &box)
        : cornerIndex(box.cornerIndex) {
    setNBoxes(box.nBoxes);
    setDerivedParameters();
}

template<int D>
BoundingBox<D> &BoundingBox<D>::operator=(const BoundingBox<D> &box) {
    if (&box != this) {
        this->cornerIndex = box.cornerIndex;
        setNBoxes(box.nBoxes);
        setDerivedParameters();
    }
    return *this;
}

template<int D>
void BoundingBox<D>::setNBoxes(const int *nb) {
    this->nBoxes[D] = 1;
    for (int d = 0; d < D; d++) {
        if (nb == 0) {
            this->nBoxes[d] = 1;
        } else {
            if (nb[d] <= 0) MSG_ERROR("Invalid box size");
            this->nBoxes[d] = nb[d];
            this->nBoxes[D] *= this->nBoxes[d];
        }
    }
}

template<int D>
void BoundingBox<D>::setDerivedParameters() {
    assert(this->nBoxes[D] > 0);
    int scale = this->cornerIndex.getScale();
    const int *l = this->cornerIndex.getTranslation();
    this->unitLength = pow(2.0, -scale);
    for (int d = 0; d < D; d++) {
        assert(this->nBoxes[d] > 0);
        this->boxLengths[d] = this->unitLength * this->nBoxes[d];
        this->lowerBounds[d] = l[d] * this->unitLength;
        this->upperBounds[d] = this->lowerBounds[d] + this->boxLengths[d];
    }
}

template<int D>
NodeIndex<D> BoundingBox<D>::getNodeIndex(const double *r) const {
    assert(r != 0);
    int idx[D];
    for (int d = 0; d < D; d++) {
        double x = r[d];
        assert(x >= this->lowerBounds[d]);
        assert(x < this->upperBounds[d]);
        double div = (x - this->lowerBounds[d]) / this->unitLength;
        double iint;
        modf(div,&iint);
        idx[d] = (int) iint;
    }

    const int *cl = this->cornerIndex.getTranslation();

    int l[D];
    for (int d = 0; d < D; d++) {
        l[d] = idx[d] + cl[d];
    }

    int n = getScale();
    NodeIndex<D> nIdx(n, l);
    return nIdx;
}

template<>
NodeIndex<1> BoundingBox<1>::getNodeIndex(int bIdx) const {
    assert(bIdx >= 0 and bIdx <= nBoxes[1]);
    int n = getScale();
    int cl = this->cornerIndex.getTranslation(0);
    int l = bIdx + cl;
    NodeIndex<1> nIdx(n, &l);
    return nIdx;
}

template<int D>
NodeIndex<D> BoundingBox<D>::getNodeIndex(int bIdx) const {
    assert(bIdx >= 0 and bIdx <= nBoxes[D]);
    int l[D];
    for (int d = D - 1; d >= 0; d--) {
        int ncells = 1;
        for (int i = 0; i < d; i++) {
            ncells *= this->nBoxes[i];
        }
        double div = bIdx / ncells;
        double iint;
        modf(div, &iint);
        l[d] = (int) iint;
        bIdx -= ncells * l[d];
    }

    int n = getScale();
    const int *cl = this->cornerIndex.getTranslation();
    for (int d = 0; d < D; d++) {
        l[d] += cl[d];
    }
    NodeIndex<D> nIdx(n, l);
    return nIdx;
}

template<>
int BoundingBox<1>::getBoxIndex(const double *r) const {
    assert(r != 0);
    double x = r[0];
    if (x < this->lowerBounds[0]) return -1;
    if (x >= this->upperBounds[0]) return -1;
    double div = (x - this->lowerBounds[0]) / this->unitLength;
    double iint;
    modf(div,&iint);
    return (int) iint;
}

template<int D>
int BoundingBox<D>::getBoxIndex(const double *r) const {
    assert(r != 0);
    int idx[D];
    for (int d = 0; d < D; d++) {
        double x = r[d];
        if (x < this->lowerBounds[d]) return -1;
        if (x >= this->upperBounds[d]) return -1;
        double div = (x - this->lowerBounds[d]) / this->unitLength;
        double iint;
        modf(div,&iint);
        idx[d] = (int) iint;
    }

    int bIdx = 0;
    for (int i = D - 1; i >= 0; i--) {
        int ncells = 1;
        for (int j = 0; j < i; j++) {
            ncells *= this->nBoxes[j];
        }
        bIdx += ncells * idx[i];
    }
    return bIdx;
}

template<>
int BoundingBox<1>::getBoxIndex(const NodeIndex<1> &nIdx) const {
    int n = nIdx.getScale();
    int l = nIdx.getTranslation(0);
    int cn = this->cornerIndex.getScale();
    int cl = this->cornerIndex.getTranslation(0);
    int relScale = n - cn;
    if (relScale >= 0) {
        return (l >> relScale) - cl;
    } else {
        return -1;
    }
}

template<int D>
int BoundingBox<D>::getBoxIndex(const NodeIndex<D> &nIdx) const {
    int n = nIdx.getScale();
    int cn = this->cornerIndex.getScale();
    const int *l = nIdx.getTranslation();
    const int *cl = this->cornerIndex.getTranslation();
    int relScale = n - cn;
    if (relScale < 0) return -1;

    int bIdx = 0;
    for (int d = D - 1; d >= 0; d--) {
        int ncells = 1;
        for (int i = 0; i < d; i++) {
            ncells *= this->nBoxes[i];
        }
        int reqTransl = (l[d] >> relScale) - cl[d];
        bIdx += ncells * reqTransl;
    }
    return bIdx;
}

template class BoundingBox<1>;
template class BoundingBox<2>;
template class BoundingBox<3>;
