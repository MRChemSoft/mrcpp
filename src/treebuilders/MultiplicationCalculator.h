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
 * @file
 * @brief Node-wise calculator for the pointwise product of multiple trees.
 *
 * @details
 * This header defines @ref mrcpp::MultiplicationCalculator, a
 * @ref mrcpp::TreeCalculator that builds the coefficients of an output node as
 * the **product** of the corresponding nodes from a vector of input
 * @ref mrcpp::FunctionTree objects.
 *
 * For each output node:
 * 1. The node index is taken from the output tree.
 * 2. For every input tree, the corresponding node is **materialized** (created
 *    if absent) via `getNode(idx)`, reconstructed to scaling space
 *    (`mwTransform(Reconstruction)`) and brought to coefficient vector form
 *    (`cvTransform(Forward)`).
 * 3. Coefficients are multiplied element-wise across all inputs, optionally
 *    applying complex conjugation according to the calculator's configuration
 *    and per-tree flags.
 * 4. The output node is transformed back (`cvTransform(Backward)`), compressed
 *    (`mwTransform(Compression)`), flagged as having coefficients, and its
 *    norms are updated.
 *
 * The calculator supports both real and complex coefficient types. When
 * `T == ComplexDouble`, conjugation is applied if either:
 * - the specific input tree advertises `conjugate()==true`, or
 * - the calculator is constructed with `conjugate=true` and the input is the
 *   **first** factor (used to implement bra–ket style products).
 */

#include "TreeCalculator.h"
#include "trees/FunctionTreeVector.h"

namespace mrcpp {

/**
 * @class MultiplicationCalculator
 * @brief Computes the pointwise product of several input trees into an output tree.
 *
 * @tparam D Spatial dimension (e.g., 1, 2, 3).
 * @tparam T Coefficient scalar type (`double` or `ComplexDouble`).
 *
 * @details
 * Let \f$\{f_i\}\f$ be the input trees (with optional scalar prefactors
 * provided externally via `get_coef(prod_vec, i)`) and let
 * \f$g = \prod_i f_i\f$ denote the pointwise product. This calculator fills
 * the coefficients of the output node corresponding to a given
 * @ref NodeIndex by:
 *
 * \f[
 *   \mathbf{c}^{(g)} \leftarrow
 *   \prod_i \left( c_i \; \mathbf{c}^{(f_i)} \right),
 * \f]
 *
 * where \f$c_i\f$ is the scalar multiplier returned by `get_coef` and
 * \f$\mathbf{c}^{(f_i)}\f$ are the (reconstructed, forward-transformed)
 * coefficients of the input node. When `T` is complex, \f$\mathbf{c}^{(f_i)}\f$
 * may be conjugated as described below.
 *
 * ### Conjugation rules (complex case)
 * - If `func_i.conjugate()` is true, that input’s coefficients are conjugated.
 * - Additionally, if the calculator is constructed with `conjugate=true`,
 *   the **first** input (index 0) is conjugated. The two conditions are XOR’d
 *   (`xor`), so a per-tree conjugation flag can cancel the global one.
 *
 * @note Missing input nodes are generated on demand by `FunctionTree::getNode`.
 */
template <int D, typename T>
class MultiplicationCalculator final : public TreeCalculator<D, T> {
public:
    /**
     * @brief Construct a product calculator.
     *
     * @param inp       Vector of input trees to be multiplied.
     * @param conjugate If `true`, apply complex conjugation to the **first**
     *                  input factor (useful for ⟨bra|ket⟩-like operations).
     *                  Ignored for real types.
     */
    MultiplicationCalculator(const FunctionTreeVector<D, T> &inp, bool conjugate = false)
            : prod_vec(inp)
            , conj(conjugate) {}

private:
    /// Collection of input trees and (optionally) their scalar prefactors.
    FunctionTreeVector<D, T> prod_vec;

    /// Global conjugation switch for the first factor (complex types only).
    bool conj;

    /**
     * @brief Compute coefficients for one output node by multiplying inputs.
     *
     * @param node_o Output node to be filled/updated.
     *
     * @details
     * Steps performed:
     * 1. Initialize output coefficients to unity.
     * 2. For each input tree `i`:
     *    - Fetch scalar factor `c_i = get_coef(prod_vec, i)`.
     *    - Materialize copy of matching input node (`getNode(idx)`),
     *      reconstruct (`mwTransform(Reconstruction)`),
     *      and forward transform to coefficient space (`cvTransform(Forward)`).
     *    - Multiply output coefficients element-wise by
     *      `c_i * coefs_i[j]` (or `c_i * conj(coefs_i[j])` per rules above).
     * 3. Transform output node back (`cvTransform(Backward)`),
     *    compress (`mwTransform(Compression)`), mark as having coefficients,
     *    and update norms.
     *
     * @note Uses helper functions `get_func(prod_vec, i)` and
     *       `get_coef(prod_vec, i)` provided by the @ref FunctionTreeVector API.
     */
    void calcNode(MWNode<D, T> &node_o) {
        const NodeIndex<D> &idx = node_o.getNodeIndex();
        T *coefs_o = node_o.getCoefs();

        // 1) Initialize output coefficients to multiplicative identity.
        for (int j = 0; j < node_o.getNCoefs(); j++) { coefs_o[j] = static_cast<T>(1.0); }

        // 2) Multiply contributions from each input factor.
        for (int i = 0; i < this->prod_vec.size(); i++) {
            T c_i = get_coef(this->prod_vec, i);
            FunctionTree<D, T> &func_i = get_func(this->prod_vec, i);

            // Materialize and prepare input node coefficients.
            MWNode<D, T> node_i = func_i.getNode(idx); // copy/materialize
            node_i.mwTransform(Reconstruction);
            node_i.cvTransform(Forward);

            const T *coefs_i = node_i.getCoefs();
            int n_coefs = node_i.getNCoefs();

            if constexpr (std::is_same<T, ComplexDouble>::value) {
                // Conjugate rule: per-tree flag XOR global-first-factor flag.
                bool do_conj = func_i.conjugate() xor (conj && i == 0);
                if (do_conj) {
                    for (int j = 0; j < n_coefs; j++) { coefs_o[j] *= c_i * std::conj(coefs_i[j]); }
                } else {
                    for (int j = 0; j < n_coefs; j++) { coefs_o[j] *= c_i * coefs_i[j]; }
                }
            } else {
                for (int j = 0; j < n_coefs; j++) { coefs_o[j] *= c_i * coefs_i[j]; }
            }
        }

        // 3) Finalize output node state.
        node_o.cvTransform(Backward);
        node_o.mwTransform(Compression);
        node_o.setHasCoefs();
        node_o.calcNorms();
    }
};

} // namespace mrcpp