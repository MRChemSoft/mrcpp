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

#include "ConvolutionOperator.h"

#include "core/InterpolatingBasis.h"
#include "core/LegendreBasis.h"

#include "functions/GaussExp.h"
#include "functions/Gaussian.h"

#include "treebuilders/CrossCorrelationCalculator.h"
#include "treebuilders/OperatorAdaptor.h"
#include "treebuilders/TreeBuilder.h"
#include "treebuilders/grid.h"
#include "treebuilders/project.h"

#include "trees/BandWidth.h"
#include "trees/FunctionTreeVector.h"
#include "trees/OperatorTree.h"

#include "utils/Printer.h"
#include "utils/Timer.h"
#include "utils/math_utils.h"

namespace mrcpp {

/**
 * @brief Construct a separable D-dimensional convolution operator from a 1D Gaussian expansion.
 *
 * The input kernel is a 1D Gaussian expansion (sum of Gauss terms). The implementation
 * projects each 1D Gaussian to a 1D function tree and then uses cross-correlations to
 * lift it into a 2D operator block; the full D-dimensional operator is assembled as a
 * separable product of these 1D blocks. The final separable rank equals kernel.size().
 *
 * @tparam D               Spatial dimension of the target operator.
 * @param mra              Multiresolution analysis defining the D-dimensional domain/basis.
 * @param kernel           1D Gaussian expansion whose terms become the separable factors.
 * @param prec             Target build precision for the operator (used for both kernel
 *                         projection and operator assembly with a small safety split).
 *
 * @details
 * Internally we choose `k_prec = prec / 10` (stricter) for fitting each 1D kernel term,
 * and `o_prec = prec` for assembling/operatorization, to keep the composed error within
 * the requested tolerance. After all factors are built, `initOperExp(kernel.size())`
 * finalizes the separable structure.
 */
template <int D>
ConvolutionOperator<D>::ConvolutionOperator(const MultiResolutionAnalysis<D> &mra,
                                            GaussExp<1> &kernel,
                                            double prec)
        : MWOperator<D>(mra, mra.getRootScale(), -10) {
    int oldlevel = Printer::setPrintLevel(0);

    this->setBuildPrec(prec);
    auto o_prec = prec;           // precision for operator assembly (2D blocks, transforms)
    auto k_prec = prec / 10.0;    // stricter precision for 1D kernel projection
    initialize(kernel, k_prec, o_prec);
    this->initOperExp(kernel.size()); // separable rank = number of kernel terms

    Printer::setPrintLevel(oldlevel);
}

/**
 * @brief Construct a convolution operator with explicit root scale and reach.
 *
 * This variant allows overriding the default operator root scale and reach (stencil
 * half-width in levels). The rest of the pipeline is identical to the other ctor:
 * build 1D kernel function trees, lift to 2D operator blocks by cross-correlation,
 * transform/collapse, then finalize the separable expansion.
 *
 * @param mra      D-dimensional MRA.
 * @param kernel   1D Gaussian expansion (rank = kernel.size()).
 * @param prec     Target build precision; we use `k_prec = prec / 100` here to be extra
 *                 conservative when the reach is user-controlled.
 * @param root     Operator root scale.
 * @param reach    Operator reach (levels outward from root). Negative = auto from box.
 */
template <int D>
ConvolutionOperator<D>::ConvolutionOperator(const MultiResolutionAnalysis<D> &mra,
                                            GaussExp<1> &kernel,
                                            double prec,
                                            int root,
                                            int reach)
        : MWOperator<D>(mra, root, reach) {
    int oldlevel = Printer::setPrintLevel(0);

    this->setBuildPrec(prec);
    auto o_prec = prec;
    auto k_prec = prec / 100.0;  // even tighter kernel fit when reach is custom
    initialize(kernel, k_prec, o_prec);
    this->initOperExp(kernel.size());

    Printer::setPrintLevel(oldlevel);
}

/**
 * @brief Core build routine: from 1D Gaussian terms → 1D function trees → 2D operator blocks.
 *
 * Steps per Gaussian term:
 *  1) **Rescaling for D-dimensional separability**: adjust the coefficient so that the product
 *     of D identical 1D factors yields the original 1D amplitude in D-D composition.
 *     Concretely: `coef ← sign(coef) * |coef|^{1/D}`.
 *  2) **Projection to a 1D function tree** (@ref FunctionTree): build an empty grid sized for
 *     narrow Gaussians (build_grid), then project the analytic Gaussian into the tree with
 *     requested kernel precision `k_prec` (project).
 *  3) **Lifting to a 2D operator**: create a @ref CrossCorrelationCalculator from the 1D tree,
 *     then use @ref TreeBuilder to expand a 2D operator tree through cross-correlations
 *     (effectively computing the correlation between basis functions along one axis).
 *  4) **Wavelet transform & caching**: bottom-up transform, compute norms, and set up node cache.
 *
 * The produced 2D blocks are stored into `raw_exp`. Higher-level code composes D-D separable
 * operators from these blocks (e.g., via @ref MWOperator’s machinery).
 *
 * @param kernel 1D Gaussian expansion (input rank).
 * @param k_prec Precision for kernel projection to 1D trees.
 * @param o_prec Precision for operator building / assembly.
 */
template <int D>
void ConvolutionOperator<D>::initialize(GaussExp<1> &kernel, double k_prec, double o_prec) {
    // Build the auxiliary 1D MRA for the kernel and fetch the D-D operator MRA
    auto k_mra = this->getKernelMRA();
    auto o_mra = this->getOperatorMRA();

    TreeBuilder<2> builder;                       // builds 2D operator trees from calculators
    OperatorAdaptor adaptor(o_prec, o_mra.getMaxScale()); // controls assembly precision / scale cap

    for (int i = 0; i < kernel.size(); i++) {
        // --- (1) Adjust coefficient for separable D-fold composition ---
        auto *k_func = kernel.getFunc(i).copy();
        // Raise absolute coefficient to 1/D and reapply sign to preserve signed kernels
        k_func->setCoef(std::copysign(std::pow(std::abs(k_func->getCoef()), 1.0 / D),
                                      k_func->getCoef()));

        // --- (2) Project analytic Gaussian to a 1D function tree ---
        FunctionTree<1> k_tree(k_mra);
        mrcpp::build_grid(k_tree, *k_func);      // Prepare an empty grid (fine where Gaussian is narrow)
        mrcpp::project(k_prec, k_tree, *k_func); // Fit the Gaussian into the 1D multiresolution basis
        delete k_func;                            // No longer needed; k_tree holds the discretization

        // --- (3) Lift to a 2D operator via cross-correlation ---
        CrossCorrelationCalculator calculator(k_tree);
        auto o_tree = std::make_unique<OperatorTree>(o_mra, o_prec);
        builder.build(*o_tree, calculator, adaptor, -1); // Dense 2D operator block in MW format

        // --- (4) Transform, normalize, and cache for application ---
        Timer trans_t;
        o_tree->mwTransform(BottomUp);  // move to MW (scaling+wavelet) representation efficiently
        o_tree->calcSquareNorm();       // useful for diagnostics / thresholding
        o_tree->setupOperNodeCache();   // enable fast repeated applications
        print::time(10, "Time transform", trans_t);
        print::separator(10, ' ');

        this->raw_exp.push_back(std::move(o_tree));
    }
}

/**
 * @brief Build a 1D MRA used to discretize the kernel factors.
 *
 * The kernel MRA mirrors the scaling family used by the D-D operator MRA:
 *  - If the operator uses an interpolating basis of order s, the kernel basis is
 *    chosen as InterpolatingBasis with order `2*s + 1`.
 *  - If Legendre, we similarly pick a LegendreBasis of order `2*s + 1`.
 *
 * The box extent (reach) is derived from the D-D world box unless an explicit
 * operator reach was set. The same uniform scaling factor is used.
 *
 * @return A standalone 1D @ref MultiResolutionAnalysis matching the operator’s scaling family.
 */
template <int D>
MultiResolutionAnalysis<1> ConvolutionOperator<D>::getKernelMRA() const {
    const BoundingBox<D> &box = this->MRA.getWorldBox();
    const ScalingBasis &basis = this->MRA.getScalingBasis();

    // Choose a kernel basis compatible with the operator basis.
    int type = basis.getScalingType();
    int kern_order = 2 * basis.getScalingOrder() + 1; // (2s+1) ensures adequate quadrature/correlation support

    ScalingBasis *kern_basis = nullptr;
    if (type == Interpol) {
        kern_basis = new InterpolatingBasis(kern_order);
    } else if (type == Legendre) {
        kern_basis = new LegendreBasis(kern_order);
    } else {
        MSG_ABORT("Invalid scaling type");
    }

    // Kernel root = operator root; reach defaults to the maximum box extent if negative.
    int root = this->oper_root;
    int reach = this->oper_reach + 1; // +1 because the 1D kernel must cover neighbors used by correlations
    if (reach < 0) {
        for (int i = 0; i < D; i++) {
            if (box.size(i) > reach) reach = box.size(i);
        }
    }

    // Build a 1D bounding box centered at zero:
    // levels from -reach to +reach (total 2*reach) at the operator root scale.
    auto start_l = std::array<int, 1>{-reach};
    auto tot_l   = std::array<int, 1>{2 * reach};
    // Uniform scaling factor (operators are implemented for uniform scales only)
    auto sf = std::array<double, 1>{box.getScalingFactor(0)};
    BoundingBox<1> kern_box(root, start_l, tot_l, sf);

    MultiResolutionAnalysis<1> kern_mra(kern_box, *kern_basis);
    delete kern_basis;
    return kern_mra;
}

// Explicit template instantiations for the supported dimensionalities.
template class ConvolutionOperator<1>;
template class ConvolutionOperator<2>;
template class ConvolutionOperator<3>;

} // namespace mrcpp
