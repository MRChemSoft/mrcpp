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

/*
 * \breif Simple storage class for scale and translation indexes.
 * The usefulness of the class becomes evident when examining
 * the parallel algorithms for projection & friends.
 */

#pragma once

#include <iostream>
#include <iomanip>

namespace mrcpp {

template <int D> class NodeIndex final {
public:
    // regular constructors
    NodeIndex(int n = 0, const std::array<int, D> &l = {})
            : N(static_cast<short int>(n))
            , L(l) {}

    // relative constructors
    NodeIndex<D> parent() const {
        std::array<int, D> l;
        for (int d = 0; d < D; d++) l[d] = (this->L[d] < 0) ? (this->L[d] - 1) / 2 : this->L[d] / 2;
        return NodeIndex<D>(this->N - 1, l);
    }
    NodeIndex<D> child(int cIdx) const {
        std::array<int, D> l;
        for (int d = 0; d < D; d++) l[d] = (2 * this->L[d]) + ((cIdx >> d) & 1);
        return NodeIndex<D>(this->N + 1, l);
    }

    // comparisons
    bool operator!=(const NodeIndex<D> &idx) const { return not(*this == idx); }
    bool operator==(const NodeIndex<D> &idx) const {
        bool out = (this->N == idx.N);
        for (int d = 0; d < D; d++) out &= (this->L[d] == idx.L[d]);
        return out;
    }
    // defines an order of the nodes (allows to use std::map)
    bool operator<(const NodeIndex<D> &idy) const {
        const NodeIndex<D> &idx = *this;
        if (idx.N != idy.N) return idx.N < idy.N;
        if (idx.L[0] != idy.L[0] or D < 2) return idx.L[0] < idy.L[0];
        if (idx.L[1] != idy.L[1] or D < 3) return idx.L[1] < idy.L[1];
        return idx.L[2] < idy.L[2];
    }

    // setters
    void setScale(int n) { this->N = static_cast<short int>(n); }
    void setTranslation(const std::array<int, D> &l) { this->L = l; }

    // value getters
    int getScale() const { return this->N; }
    int getTranslation(int d) const { return this->L[d]; }
    std::array<int, D> getTranslation() const { return this->L; }

    // reference getters
    int &operator[](int d) { return this->L[d]; }
    const int &operator[](int d) const { return this->L[d]; }

    std::ostream &print(std::ostream &o) const {
        o << "[ " << std::setw(3) << this->N << " | ";
        for (int d = 0; d < D - 1; d++) o << std::setw(4) << this->L[d] << ", ";
        o << std::setw(4) << this->L[D - 1] << "]";
        return o;
    }

private:
    short int N{0};          ///< Length scale index 2^N
    std::array<int, D> L{}; ///< Translation index [x,y,z,...]
};

/** @brief ostream printer */
template <int D> std::ostream &operator<<(std::ostream &o, const NodeIndex<D> &idx) {
    return idx.print(o);
}

/** @brief Check whether indices are directly related (not sibling) */
template <int D> bool related(const NodeIndex<D> &a, const NodeIndex<D> &b) {
    const auto &sr = (a.getScale() < b.getScale()) ? a : b;
    const auto &jr = (a.getScale() >= b.getScale()) ? a : b;
    auto rel_scale = jr.getScale() - sr.getScale();

    bool related = true;
    for (int d = 0; d < D; d++) related &= (sr[d] == (jr[d] >> rel_scale));
    return related;
}

/** @brief Check whether indices are siblings, i.e. same parent */
template <int D> bool siblings(const NodeIndex<D> &a, const NodeIndex<D> &b) {
    return (a.parent() == b.parent());
}

} // namespace mrcpp
