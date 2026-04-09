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
 * @brief Node-wise power transform calculator for multiresolution trees.
 *
 * @details
 * This header defines @ref mrcpp::PowerCalculator, a concrete
 * @ref TreeCalculator that raises the coefficients of an input function tree
 * to a scalar power, writing results into an output tree during the standard
 * calculator traversal.
 *
 * The transform is applied **locally per node** in scaling space:
 *  1. Fetch/copy the corresponding input node (creating it if missing).
 *  2. Apply multiresolution reconstruction (wavelet → scaling) and then a
 *     forward coefficient transform to obtain coefficient values suitable
 *     for pointwise operations.
 *  3. Compute @c coefs_out[j] = pow(coefs_in[j], power) for each coefficient.
 *  4. Apply the inverse coefficient transform and multiresolution compression.
 *  5. Mark coefficients present and update node norms.
 *
 * The scalar @ref power is a real number. For complex-valued trees,
 * @c std::pow(std::complex<>, double) is used.
 */

#include "TreeCalculator.h"

namespace mrcpp {

/**
 * @class PowerCalculator
 * @brief Raises node coefficients of an input tree to a fixed power.
 *
 * @tparam D Spatial dimension of the tree (e.g., 1, 2, or 3).
 * @tparam T Coefficient scalar type (e.g., @c double or @c ComplexDouble).
 *
 * @details
 * This calculator implements a **pointwise power** operation in coefficient
 * space. For each visited node in the output tree, it pulls the corresponding
 * node from the input tree (creating it if necessary), reconstructs to scaling
 * space, and applies:
 * @f[
 *   \forall j:\quad c^{\text{out}}_j \leftarrow \big(c^{\text{in}}_j\big)^{\,p}
 * @f]
 * where @f$p=@ref power@f$.
 *
 * ### Precision & grid handling
 * The class delegates traversal, splitting, and precision handling to the
 * base @ref TreeCalculator. It performs only the node-local algebra and the
 * required forward/backward transforms.
 *
 * ### Complex inputs
 * For complex-valued @p T, the standard overload @c std::pow(T,double) is used.
 * Note that if @p T is real and coefficients are negative while @ref power is
 * non-integer, the result may be @c NaN; this calculator does not alter that
 * behavior.
 */
template <int D, typename T>
class PowerCalculator final : public TreeCalculator<D, T> {
public:
    /**
     * @brief Construct a power calculator for a given input tree.
     *
     * @param inp  Reference to the input function tree whose node coefficients
     *             are the base values in the power operation.
     * @param pow  Exponent @f$p@f$ to apply pointwise to all node coefficients.
     *
     * @note The calculator does not own @p inp; the caller must ensure that
     *       the referenced tree remains valid for the calculator's lifetime.
     */
    PowerCalculator(FunctionTree<D, T> &inp, double pow)
            : power(pow)
            , func(&inp) {}

private:
    /**
     * @brief Scalar exponent used in @c std::pow for all coefficients.
     */
    double power;

    /**
     * @brief Non-owning pointer to the input function tree.
     */
    FunctionTree<D, T> *func;

    /**
     * @brief Node-level power application.
     *
     * @param node_o Output node whose coefficients are overwritten with
     *               @f$\big(c^{\text{in}}\big)^{\,p}@f$ at the corresponding
     *               location in the input tree.
     *
     * @details
     * Steps performed:
     *  - Retrieve the corresponding input node @c node_i from @ref func
     *    using the same @ref NodeIndex (this may create a missing node).
     *  - @c node_i.mwTransform(Reconstruction) to convert wavelet → scaling.
     *  - @c node_i.cvTransform(Forward) to access coefficient array in the
     *    appropriate local basis.
     *  - For each coefficient index @c j, compute:
     *    @code
     *    coefs_o[j] = std::pow(coefs_i[j], power);
     *    @endcode
     *  - Apply @c node_o.cvTransform(Backward) and
     *    @c node_o.mwTransform(Compression) to restore representation.
     *  - Set the "has coefficients" flag and refresh node norms.
     *
     * @warning If @p T is a real type and @ref power is non-integer, negative
     *          input coefficients can lead to @c NaN. This is the standard
     *          behavior of @c std::pow and is not intercepted here.
     */
    void calcNode(MWNode<D, T> &node_o) override {
        const NodeIndex<D> &idx = node_o.getNodeIndex();
        int n_coefs = node_o.getNCoefs();
        T *coefs_o = node_o.getCoefs();

        // Generate/copy input node at the same index.
        MWNode<D, T> node_i = func->getNode(idx);
        node_i.mwTransform(Reconstruction);
        node_i.cvTransform(Forward);

        const T *coefs_i = node_i.getCoefs();
        for (int j = 0; j < n_coefs; j++) { coefs_o[j] = std::pow(coefs_i[j], this->power); }

        node_o.cvTransform(Backward);
        node_o.mwTransform(Compression);
        node_o.setHasCoefs();
        node_o.calcNorms();
    }
};

} // namespace mrcpp