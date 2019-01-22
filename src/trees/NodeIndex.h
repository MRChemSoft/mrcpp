/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2019 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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

/*
 * \breif Simple storage class for scale and translation indexes.
 * The usefulness of the class becomes evident when examining
 * the parallel algorithms for projection & friends.
 */

#pragma once

#include <iostream>

#include "mrcpp_declarations.h"

namespace mrcpp {

template<int D>
class NodeIndex final {
public:
    NodeIndex(int n = 0, const int *l = nullptr);
    NodeIndex(const NodeIndex<D> &idx);
    NodeIndex(const NodeIndex<D> &pIdx, int cIdx);

    inline NodeIndex<D>& operator=(const NodeIndex<D> &idx);
    inline bool operator==(const NodeIndex<D> &idx) const;
    inline bool operator!=(const NodeIndex<D> &idx) const;

    void setScale(int n) { this->N = (short int) n; }
    inline void setTranslation(const int *l);

    int getScale() const { return this->N; }
    int getTranslation(int d) const { return this->L[d]; }
    int *getTranslation() { return this->L; }
    const int *getTranslation() const { return this->L; }

    friend std::ostream& operator<<(std::ostream &o, const NodeIndex<D> &idx) { return idx.print(o); }
    friend class NodeIndexComp<D>;

private:
    short int N;
    int L[D];

    std::ostream& print(std::ostream &o) const;
};

template<int D>
NodeIndex<D>::NodeIndex(int n, const int *l) {
    this->N = (short int) n;
    setTranslation(l);
}

template<int D>
NodeIndex<D>::NodeIndex(const NodeIndex<D> &idx) {
    this->N = idx.N;
    setTranslation(idx.L);
}

template<int D>
NodeIndex<D>::NodeIndex(const NodeIndex<D> &pIdx, int cIdx) {
    this->N = pIdx.N + 1;
    const int *l = pIdx.getTranslation();
    for (int d = 0; d < D; d++) {
        this->L[d] = (2 * l[d]) + ((cIdx >> d) & 1);
    }
}

template<int D>
NodeIndex<D>& NodeIndex<D>::operator=(const NodeIndex<D> &idx) {
    if (&idx == this) {
        return *this;
    }
    this->N = idx.N;
    setTranslation(idx.L);
    return *this;
}

template<int D>
void NodeIndex<D>::setTranslation(const int *l) {
    for (int d = 0; d < D; d++) {
        if (l != nullptr) {
            this->L[d] = l[d];
        } else {
            this->L[d] = 0;
        }
    }
}

template<int D>
bool NodeIndex<D>::operator==(const NodeIndex<D> &idx) const {
    if (this->N != idx.N) return false;
    for (int d = 0; d < D; d++) {
        if (this->L[d] != idx.L[d]) return false;
    }
    return true;
}

template<int D>
bool NodeIndex<D>::operator!=(const NodeIndex<D> &idx) const {
    if (this->N != idx.N) return true;
    for (int d = 0; d < D; d++) {
        if (this->L[d] != idx.L[d]) return true;
    }
    return false;
}

template<int D>
std::ostream& NodeIndex<D>::print(std::ostream &o) const {
    o << "[ " << this->N << " | ";
    for (int d = 0; d < D - 1; d++) {
        o << this->L[d] << ", ";
    }
    o << this->L[D - 1] << "]";
    return o;
}


template<int D>
class NodeIndexComp {
public:
    bool operator()(const NodeIndex<D> &a, const NodeIndex<D> &b) const {
        if (a.N < b.N) {
            return true;
        }
        if (a.N > b.N) {
            return false;
        }
        for (int d = 0; d < D; d++) {
            if (a.L[d] == b.L[d]) {
                continue;
            }
            if (a.L[d] < b.L[d]) {
                return true;
            }
            return false;
        }
        return false;
    }

    bool operator()(const NodeIndex<D> *a, const NodeIndex<D> *b) const {
        if (a->N < b->N) {
            return true;
        }
        if (a->N > b->N) {
            return false;
        }
        for (int d = 0; d < D; d++) {
            if (a->L[d] == b->L[d]) {
                continue;
            }
            if (a->L[d] < b->L[d]) {
                return true;
            }
            return false;
        }
        return false;
    }
};

}
