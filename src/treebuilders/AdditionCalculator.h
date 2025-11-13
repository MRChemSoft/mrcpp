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
/**
 * @file AdditionCalculator.h
 * @brief Node-wise accumulator used during adaptive construction to sum
 *        multiresolution (MW) functions with optional conjugation.
 *
 * @details
 * This header defines #mrcpp::AdditionCalculator, a lightweight
 * #mrcpp::TreeCalculator that, for each target node, fetches the
 * corresponding node from every input function in a
 * #mrcpp::FunctionTreeVector and accumulates a weighted sum of their
 * coefficients. No refinement policy is implemented here; pair this
 * calculator with a #mrcpp::TreeBuilder and a suitable #mrcpp::TreeAdaptor.
 *
 * Complex handling:
 * - For complex `T`, each term uses either the raw coefficients or their
 *   complex conjugate according to the XOR of the input tree's own
 *   `conjugate()` flag and the calculator-wide `conj` flag.
 */

#include <type_traits>   // std::is_same
#include <complex>       // std::conj

#include "TreeCalculator.h"
#include "trees/FunctionTreeVector.h"

namespace mrcpp {

/**
 * @class AdditionCalculator
 * @brief Node-wise accumulator for adaptive sums of multiresolution functions.
 *
 * @tparam D Spatial dimension (1, 2, or 3).
 * @tparam T Coefficient type (e.g., `double`, `ComplexDouble`).
 *
 * @details
 * For each target node \p node_o (identified by its NodeIndex), this calculator
 * gathers the corresponding node from every input function in a
 * #mrcpp::FunctionTreeVector and accumulates the weighted coefficients into
 * \p node_o:
 *
 * \f[
 *   \mathbf{c}_o \;=\; \sum_i \alpha_i \,\mathbf{c}_i .
 * \f]
 *
 * If \p T is complex, optional conjugation is applied according to the XOR of
 * the per-tree conjugation flag and the calculator-wide @ref conj flag; i.e.
 * a term uses \f$\overline{\mathbf{c}_i}\f$ iff exactly one of the two flags is set.
 *
 * This class performs **no grid refinement** or transforms; it only writes
 * coefficients, marks presence, and updates node norms. It is intended to be
 * used inside the adaptive loop driven by #mrcpp::TreeBuilder together with an
 * appropriate adaptor.
 */
template <int D, typename T>
class AdditionCalculator final : public TreeCalculator<D, T> {
public:
    /**
     * @brief Construct an addition calculator over a set of input trees.
     *
     * @param[in] inp       Vector of (coefficient, tree) pairs to be summed.
     * @param[in] conjugate Global conjugation toggle for complex types. For
     *                      each input tree, the effective conjugation applied
     *                      is `tree.conjugate() XOR conjugate`.
     *
     * @note All input trees are assumed to share an MRA compatible with the
     *       output tree provided to the builder.
     */
    AdditionCalculator(const FunctionTreeVector<D, T> &inp, bool conjugate = false)
            : sum_vec(inp)
            , conj(conjugate) {}

private:
    /// Vector of weighted input trees to sum.
    FunctionTreeVector<D, T> sum_vec;
    /// Global conjugation toggle for complex accumulation (see ctor docs).
    bool conj;

    /**
     * @brief Accumulate coefficients for a single output node.
     *
     * @param[in,out] node_o Target node whose coefficients are overwritten
     *                       by the weighted sum of matching input nodes.
     *
     * @details
     * Steps:
     *  1. Zero \p node_o coefficients.
     *  2. For each entry \f$(\alpha_i, f_i)\f$ in @ref sum_vec:
     *     - Fetch (and create if needed) the input node with the same index as \p node_o.
     *     - Accumulate \f$\alpha_i \cdot \mathbf{c}_i\f$ (or its conjugate for complex,
     *       following the XOR rule) into \p node_o.
     *  3. Mark coefficients present and update norms.
     *
     * No transforms are performed here; coefficients are assumed to be in the
     * same representation across all trees.
     */
    void calcNode(MWNode<D, T> &node_o) override {
        node_o.zeroCoefs();
        const NodeIndex<D> &idx = node_o.getNodeIndex();
        T *coefs_o = node_o.getCoefs();

        for (int i = 0; i < this->sum_vec.size(); i++) {
            T c_i = get_coef(this->sum_vec, i);
            FunctionTree<D, T> &func_i = get_func(this->sum_vec, i);

            // This generates the node if missing in func_i
            const MWNode<D, T> &node_i = func_i.getNode(idx);
            const T *coefs_i = node_i.getCoefs();
            int n_coefs = node_i.getNCoefs();

            if constexpr (std::is_same<T, ComplexDouble>::value) {
                const bool use_conj = (func_i.conjugate() xor conj);
                if (use_conj) {
                    for (int j = 0; j < n_coefs; j++) coefs_o[j] += c_i * std::conj(coefs_i[j]);
                } else {
                    for (int j = 0; j < n_coefs; j++) coefs_o[j] += c_i * coefs_i[j];
                }
            } else {
                for (int j = 0; j < n_coefs; j++) coefs_o[j] += c_i * coefs_i[j];
            }
        }

        node_o.setHasCoefs();
        node_o.calcNorms();
    }
};

} // namespace mrcpp