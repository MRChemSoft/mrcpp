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
 * @brief Element-wise squaring of function-tree coefficients (with optional complex conjugation).
 *
 * @details
 * This calculator evaluates one of the following pointwise operations on an input function
 * represented by a multiresolution @ref FunctionTree:
 *
 * - **Algebraic square**: \f$ g(\mathbf r) = f(\mathbf r)^2 \f$
 * - **Squared magnitude** (Hermitian square): \f$ g(\mathbf r) = f(\mathbf r)\,f(\mathbf r)^* = |f(\mathbf r)|^2 \f$
 *
 * The choice is controlled by the constructor's `conjugate` flag (see below). For real
 * coefficient types the two definitions coincide.
 *
 * Implementation sketch per node:
 * 1. Pull (or generate) the input node at the same index as the output node.
 * 2. Convert to scaling coefficients (multiwavelet reconstruction), then to coefficient
 *    vector space.
 * 3. Apply the element-wise operation (square or squared magnitude).
 * 4. Transform coefficients back, compress, mark as having coefficients, and refresh norms.
 */

#include "TreeCalculator.h"

namespace mrcpp {

/**
 * @class SquareCalculator
 * @brief Per-node square / squared-magnitude operator for function trees.
 *
 * @tparam D Spatial dimension of the tree.
 * @tparam T Scalar coefficient type (e.g., `double` or `ComplexDouble`).
 *
 * @details
 * Let \f$f\f$ be the input function represented by `func`. This calculator writes
 * an output tree \f$g\f$ such that, for each node and each basis coefficient:
 *
 * - If `conjugate == false`:
 *   - Real `T`: \f$ g = f^2 \f$
 *   - Complex `T`: \f$ g = f^2 \f$
 * - If `conjugate == true`:
 *   - Real `T`: \f$ g = f^2 \f$ (same as above)
 *   - Complex `T`: \f$ g = f\,\overline{f} = |f|^2 \f$
 *
 * Additionally, if the input tree is marked internally as "conjugated" (via its
 * soft-conjugation flag), the implementation respects that view such that
 * `conjugate == true` still produces \f$|f|^2\f$ and `conjugate == false` produces
 * \f$(\overline{f})^2\f$ for complex `T`. See the truth table in @ref calcNode.
 *
 * ### Transform pipeline
 * Each output node is computed by:
 * - reconstructing the corresponding input node to scaling space
 *   (`mwTransform(Reconstruction)`),
 * - converting to coefficient space (`cvTransform(Forward)`),
 * - applying the element-wise operation on the coefficient array,
 * - mapping back (`cvTransform(Backward)`, `mwTransform(Compression)`),
 * - finalizing (`setHasCoefs()`, `calcNorms()`).
 *
 * @note The calculator is stateless across nodes and can be scheduled in parallel by
 *       the tree execution engine as long as nodes are independent.
 */
template <int D, typename T>
class SquareCalculator final : public TreeCalculator<D, T> {
public:
    /**
     * @brief Construct a square (or squared-magnitude) calculator.
     *
     * @param[in] inp        Input function tree \f$f\f$.
     * @param[in] conjugate  If `true` and `T` is complex, compute the squared magnitude
     *                       \f$|f|^2\f$ (i.e., multiply by the complex conjugate). If
     *                       `false`, compute the algebraic square \f$f^2\f$.
     *
     * @note For real `T`, `conjugate` has no effect; \f$|f|^2 = f^2\f$.
     */
    SquareCalculator(FunctionTree<D, T> &inp, bool conjugate = false)
            : func(&inp)
            , conj(conjugate) {}

private:
    /// Pointer to the input function tree \f$f\f$.
    FunctionTree<D, T> *func;
    /// Operation switch: `false` ⇒ \f$f^2\f$; `true` ⇒ (for complex) \f$|f|^2\f$.
    bool conj;

    /**
     * @brief Compute one output node by squaring the corresponding input node.
     *
     * @param[in,out] node_o Output node to be written at the current index.
     *
     * @details
     * Steps:
     * 1. Acquire a copy of the input node at the same index: `node_i = func->getNode(idx)`.
     * 2. Transform `node_i` to coefficient space (`mwTransform(Reconstruction)`,
     *    then `cvTransform(Forward)`).
     * 3. For each coefficient \f$c_j\f$:
     *    - If `T` is real: \f$c_j \leftarrow c_j^2\f$.
     *    - If `T` is complex:
     *        - Respect the input tree's soft conjugation flag (`func->conjugate()`).
     *        - Apply the following table to compute `coefs_o[j]`:
     *
     *        | `func->conjugate()` | `conj` | result                                 |
     *        |:--------------------:|:------:|:---------------------------------------|
     *        |        false         | false  | \f$c_j \cdot c_j = c_j^2\f$            |
     *        |        false         |  true  | \f$c_j \cdot \overline{c_j} = |c_j|^2\f$ |
     *        |         true         | false  | \f$\overline{c_j}\cdot \overline{c_j} = (\overline{c_j})^2\f$ |
     *        |         true         |  true  | \f$\overline{c_j}\cdot c_j = |c_j|^2\f$ |
     *
     * 4. Map the result back (`cvTransform(Backward)`, `mwTransform(Compression)`),
     *    then finalize flags and norms.
     *
     * @complexity Linear in the number of coefficients of the node: \f$O(n_{\text{coefs}})\f$.
     */
    void calcNode(MWNode<D, T> &node_o) {
        const NodeIndex<D> &idx = node_o.getNodeIndex();
        int n_coefs = node_o.getNCoefs();
        T *coefs_o = node_o.getCoefs();

        // Acquire / materialize the input node at the same index
        MWNode<D, T> node_i = func->getNode(idx); // Copy node (may generate missing nodes)
        node_i.mwTransform(Reconstruction);
        node_i.cvTransform(Forward);

        const T *coefs_i = node_i.getCoefs();

        if constexpr (std::is_same<T, ComplexDouble>::value) {
            if (func->conjugate()) {
                if (conj) {
                    // |f|^2: conj(c) * c
                    for (int j = 0; j < n_coefs; j++) { coefs_o[j] = std::conj(coefs_i[j]) * coefs_i[j]; }
                } else {
                    // (conj f)^2: conj(c) * conj(c)
                    for (int j = 0; j < n_coefs; j++) { coefs_o[j] = std::conj(coefs_i[j]) * std::conj(coefs_i[j]); }
                }
            } else {
                if (conj) {
                    // |f|^2: c * conj(c)
                    for (int j = 0; j < n_coefs; j++) { coefs_o[j] = coefs_i[j] * std::conj(coefs_i[j]); }
                } else {
                    // f^2: c * c
                    for (int j = 0; j < n_coefs; j++) { coefs_o[j] = coefs_i[j] * coefs_i[j]; }
                }
            }
        } else {
            // Real case: f^2
            for (int j = 0; j < n_coefs; j++) { coefs_o[j] = coefs_i[j] * coefs_i[j]; }
        }

        // Map back and finalize
        node_o.cvTransform(Backward);
        node_o.mwTransform(Compression);
        node_o.setHasCoefs();
        node_o.calcNorms();
    }
};

} // namespace mrcpp