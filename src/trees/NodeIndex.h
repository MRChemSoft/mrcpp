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

#pragma once

#include <array>
#include <iomanip>
#include <iostream>

namespace mrcpp {

/**
 * @class NodeIndex
 * @tparam D Spatial dimension (1, 2, or 3)
 * @brief Storage class for scale and translation indexes
 *
 * @details
 * A NodeIndex encodes the position of a node in a multiresolution tree by:
 * - an integer **scale** N (node size proportional to 2^{-N})
 * - D-dimensional **translation** vector L og integers
 *
 * Provides helpers to obtain the parent/child indices, comparisons
 * (including a strict weak ordering for associative containers), and utilities
 * to test ancestry/sibling relations
 * 
 * The usefulness of the class becomes evident when examining
 * the parallel algorithms for projection & friends
 *
 * @note
 * The scale is stored as a short integer; the translation is a D-dimensional
 * integer vector. The translation follows the standard dyadic refinement.
 */
template <int D> class NodeIndex final {
public:
    /**
     * @brief Regular constructor for NodeIndex
     * @param[in] n Scale (defaults to zero)
     * @param[in] l Translation vector with dimension D
     *
     * @details Casts n to a short int (N) and directly assigns L as l
     */
    NodeIndex(int n = 0, const std::array<int, D> &l = {})
            : N(static_cast<short int>(n))
            , L(l) {}
    /**
     * @brief Relative constructor of the parent NodeIndex
     * @return Parent NodeIndex
     *
     * @details Parents (N = N - 1) are obtained by floor rounding L/2
     */
    NodeIndex<D> parent() const {
        std::array<int, D> l;
        for (int d = 0; d < D; d++) l[d] = (this->L[d] < 0) ? (this->L[d] - 1) / 2 : this->L[d] / 2;
        return NodeIndex<D>(this->N - 1, l);
    }
    /**
     * @brief Relative constructor of child NodeIndex
     * @param cIdx Child linear index
     * @return Child NodeIndex
     * 
     * @details Children (N = N + 1) are obtained by L = 2L + b, where @c b is given by the bits of @p cIdx
     */
    NodeIndex<D> child(int cIdx) const {
        std::array<int, D> l;
        for (int d = 0; d < D; d++) l[d] = (2 * this->L[d]) + ((cIdx >> d) & 1);
        return NodeIndex<D>(this->N + 1, l);
    }

    /**
     * @brief Defines inequality operator
     * @param[in] idx NodeIndex of the comparing node
     * @return True if N and/or L are different
     */
    bool operator!=(const NodeIndex<D> &idx) const { return not(*this == idx); }
    /**
     * @brief Defines equality operator
     * @param[in] idx NodeIndex of comparing node
     * @return True if both N and L are equal
     */
    bool operator==(const NodeIndex<D> &idx) const {
        bool out = (this->N == idx.N);
        for (int d = 0; d < D; d++) out &= (this->L[d] == idx.L[d]);
        return out;
    }
    /**
     * @brief Defines comparison operator
     * @param[in] idy NodeIndex of comparing node
     * @return True if *this is smaller than idy
     *
     * @details
     * Comparison rules (by order):
     *  1. NodeIndex with smallest N is considered smallest
     *  2. NodeIndex with the first component of L be smaller is considered smaller
     *
     * @note
     * Strict weak ordering provides strict weak ordering to enables usage in std::map
     */
    bool operator<(const NodeIndex<D> &idy) const {
        const NodeIndex<D> &idx = *this;
        if (idx.N != idy.N) return idx.N < idy.N;
        if (idx.L[0] != idy.L[0] or D < 2) return idx.L[0] < idy.L[0];
        if (idx.L[1] != idy.L[1] or D < 3) return idx.L[1] < idy.L[1];
        return idx.L[2] < idy.L[2];
    }

    /*
     * Getters and setters
     */
    int getScale() const { return this->N; }                            ///< @return Scale of node
    std::array<int, D> getTranslation() const { return this->L; }       ///< @return Full translation vector
    void setScale(int n) { this->N = static_cast<short int>(n); }       ///< @param n Scale of node
    void setTranslation(const std::array<int, D> &l) { this->L = l; }   ///< @param l Translation vector of dimension D
    
    /**
     * @brief Get a specific component of translation vector, L
     * @param[in] d Index of wanted component
     * @return Translation vector component @p d
    */
    int getTranslation(int d) const { return this->L[d]; }

    /**
     * @brief Define indexing operator of translation vector, L
     * @param[in] d Index of wanted component
     * @return Translation vector component @p d
    */
    int &operator[](int d) { return this->L[d]; }

    /**
     * @brief Const version of @ref &operator[]
     * @param[in] d Index of wanted component
     * @return Translation vector component @p d
    */
    const int &operator[](int d) const { return this->L[d]; }

    /**
     * @brief Creates output stream of NodeIndex in readable format
     * @param o Output stream
     * @return A formatted version of @o
     * 
     * @details
     * Prints NodeIndex on the form "[ N | L0, L1, ... ]"
     */
    std::ostream &print(std::ostream &o) const {
        o << "[ " << std::setw(3) << this->N << " | ";
        for (int d = 0; d < D - 1; d++) o << std::setw(4) << this->L[d] << ", ";
        o << std::setw(4) << this->L[D - 1] << "]";
        return o;
    }

private:
    short int N{0};         ///< Length scale index 2^N
    std::array<int, D> L{}; ///< Translation index [x,y,z,...]
};

/**
 * @brief Defines operator for print of a @ref NodeIndex
 * @param o Output stream
 * @param[in] idx NodeIndex of wanted node
 * @return Print stream of @ref NodeIndex
 */
template <int D> std::ostream &operator<<(std::ostream &o, const NodeIndex<D> &idx) {
    return idx.print(o);
}

/**
 * @brief Check whether two NodeIndices are directly related
 * @param[in] a First NodeIndex
 * @param[in] b Second NodeIndex
 * @return True if related
 *
 * @details @p a and @p b are related if they follow the relation rules described
 * in the @ref child() and @ref parent() constructors
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
 * @brief Check whether two NodeIndices are siblings, i.e. same parent
 * @param[in] a First NodeIndex
 * @param[in] b Second NodeIndex
 * @return True if siblings
 */
template <int D> bool siblings(const NodeIndex<D> &a, const NodeIndex<D> &b) {
    return (a.parent() == b.parent());
}

} // namespace mrcpp