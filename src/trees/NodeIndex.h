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

/**
 * @file NodeIndex.h
 * @brief Compact storage for multiresolution node indices (scale and translation).
 *
 * @details
 * A NodeIndex encodes the position of a node in a multiresolution tree by:
 * - an integer **scale** \f$N\f$ (node size \f$\propto 2^{-N}\f$), and
 * - an integer **translation** vector \f$\mathbf{L}\in\mathbb{Z}^D\f$.
 *
 * The class provides helpers to obtain the parent/child indices, comparisons
 * (including a strict weak ordering for associative containers), and utilities
 * to test ancestry/sibling relations (see free functions @ref related and
 * @ref siblings below).
 */

#pragma once

#include <array>
#include <iomanip>
#include <iostream>

namespace mrcpp {

/**
 * @class NodeIndex
 * @tparam D Spatial dimension (1, 2 or 3).
 * @brief Scale–translation pair identifying a node in a MW tree.
 *
 * @details
 * The scale is stored as a short integer; the translation is a D-dimensional
 * integer vector. The translation follows the standard dyadic refinement:
 * children are obtained by doubling each component and adding the child-bit
 * extracted from the child index.
 */
template <int D> class NodeIndex final {
public:
    /**
     * @name Constructors
     * @{
     */

    /**
     * @brief Construct from scale and translation.
     * @param n Scale \f$N\f$.
     * @param l Translation vector \f$\mathbf{L}\f$ (defaults to all zeros).
     */
    NodeIndex(int n = 0, const std::array<int, D> &l = {})
            : N(static_cast<short int>(n))
            , L(l) {}

    /**
     * @brief Index of the parent node (one level coarser).
     * @return Parent index \f$(N-1, \lfloor L/2 \rfloor)\f$ with correct rounding for negatives.
     */
    NodeIndex<D> parent() const {
        std::array<int, D> l;
        for (int d = 0; d < D; d++) l[d] = (this->L[d] < 0) ? (this->L[d] - 1) / 2 : this->L[d] / 2;
        return NodeIndex<D>(this->N - 1, l);
    }

    /**
     * @brief Index of a child node (one level finer).
     * @param cIdx Child linear index in \f$[0, 2^D)\f$; bit @c d selects the offset in dimension @c d.
     * @return Child index \f$(N+1, 2L + b)\f$ with @c b given by the bits of @p cIdx.
     */
    NodeIndex<D> child(int cIdx) const {
        std::array<int, D> l;
        for (int d = 0; d < D; d++) l[d] = (2 * this->L[d]) + ((cIdx >> d) & 1);
        return NodeIndex<D>(this->N + 1, l);
    }
    /// @}

    /**
     * @name Comparisons
     * @{
     */
    /// Inequality.
    bool operator!=(const NodeIndex<D> &idx) const { return not(*this == idx); }

    /// Equality (same scale and same translation vector).
    bool operator==(const NodeIndex<D> &idx) const {
        bool out = (this->N == idx.N);
        for (int d = 0; d < D; d++) out &= (this->L[d] == idx.L[d]);
        return out;
    }

    /**
     * @brief Strict weak ordering (by scale, then lexicographically by translation).
     * @details Enables usage as key in @c std::map / @c std::set.
     */
    bool operator<(const NodeIndex<D> &idy) const {
        const NodeIndex<D> &idx = *this;
        if (idx.N != idy.N) return idx.N < idy.N;
        if (idx.L[0] != idy.L[0] or D < 2) return idx.L[0] < idy.L[0];
        if (idx.L[1] != idy.L[1] or D < 3) return idx.L[1] < idy.L[1];
        return idx.L[2] < idy.L[2];
    }
    /// @}

    /**
     * @name Setters
     * @{
     */
    /// Set the scale.
    void setScale(int n) { this->N = static_cast<short int>(n); }

    /// Set the translation vector.
    void setTranslation(const std::array<int, D> &l) { this->L = l; }
    /// @}

    /**
     * @name Getters (values)
     * @{
     */
    /// @return The scale \f$N\f$.
    int getScale() const { return this->N; }

    /// @return Component @p d of the translation vector.
    int getTranslation(int d) const { return this->L[d]; }

    /// @return Full translation vector.
    std::array<int, D> getTranslation() const { return this->L; }
    /// @}

    /**
     * @name Getters (references)
     * @{
     */
    /// Mutable access to translation component @p d.
    int &operator[](int d) { return this->L[d]; }

    /// Const access to translation component @p d.
    const int &operator[](int d) const { return this->L[d]; }
    /// @}

    /**
     * @brief Print as "[ N | L0, L1, ... ]".
     * @param o Output stream.
     * @return The stream @p o.
     */
    std::ostream &print(std::ostream &o) const {
        o << "[ " << std::setw(3) << this->N << " | ";
        for (int d = 0; d < D - 1; d++) o << std::setw(4) << this->L[d] << ", ";
        o << std::setw(4) << this->L[D - 1] << "]";
        return o;
    }

private:
    short int N{0};         ///< Length-scale index \f$N\f$ (node size \f$\propto 2^{-N}\f$).
    std::array<int, D> L{}; ///< Translation vector \f$\mathbf{L}\f$.
};

/**
 * @brief Stream inserter for @ref NodeIndex.
 * @relates NodeIndex
 */
template <int D> std::ostream &operator<<(std::ostream &o, const NodeIndex<D> &idx) {
    return idx.print(o);
}

/**
 * @brief Test if two indices are on the same branch (ancestor/descendant relation).
 * @tparam D Dimension.
 * @param a First index.
 * @param b Second index.
 * @return @c true if the coarser index equals the finer index truncated to the coarser scale.
 *
 * @details
 * Let @c sr be the shallower (coarser) of @p a and @p b, and @c jr the deeper (finer).
 * They are related if \f$\mathbf{L}_{\text{sr}} = \lfloor \mathbf{L}_{\text{jr}} / 2^{N_{\text{jr}}-N_{\text{sr}}}\rfloor\f$.
 */
template <int D> bool related(const NodeIndex<D> &a, const NodeIndex<D> &b) {
    const auto &sr = (a.getScale() < b.getScale()) ? a : b;
    const auto &jr = (a.getScale() >= b.getScale()) ? a : b;
    auto rel_scale = jr.getScale() - sr.getScale();

    bool related = true;
    for (int d = 0; d < D; d++) related &= (sr[d] == (jr[d] >> rel_scale));
    return related;
}

/**
 * @brief Test if two indices are siblings (share the same parent).
 * @tparam D Dimension.
 * @param a First index.
 * @param b Second index.
 * @return @c true if @p a.parent() == @p b.parent().
 */
template <int D> bool siblings(const NodeIndex<D> &a, const NodeIndex<D> &b) {
    return (a.parent() == b.parent());
}

} // namespace mrcpp