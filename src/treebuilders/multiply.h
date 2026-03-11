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
 * @brief High-level algebra on multiresolution function trees.
 *
 * @details
 * This header declares scalar and field operations on
 * @ref mrcpp::FunctionTree objects:
 * - continuous inner products (dot products),
 * - pointwise products of two or many trees,
 * - powers/squares (element-wise).
 *
 * Unless stated otherwise, functions honor a target accuracy parameter
 * `prec` (see each overload). Implementations typically refine/coarsen
 * grids adaptively using wavelet norms and multiresolution estimates until
 * the requested tolerance is met (or a `maxIter` cap is reached).
 *
 * ### Precision semantics
 * - `prec` is interpreted as a **relative** tolerance in an L2-like sense
 *   by default; set `absPrec=true` to treat it as an **absolute** tolerance.
 * - `maxIter < 0` means “iterate as needed”; otherwise it limits refinement
 *   passes (the function may exit early with a looser error).
 *
 * ### Conjugation semantics (complex trees)
 * When `T` is complex, some overloads accept `conjugate=true` to apply
 * complex conjugation to the first factor (bra–ket convention), yielding
 * products like \f$f \cdot \overline{g}\f$ or \f{|f|^2}\f.
 */

#include "trees/FunctionTreeVector.h"

namespace mrcpp {

template <int D, typename T> class RepresentableFunction;
template <int D, typename T> class FunctionTree;

/**
 * @brief Continuous inner product \f$\langle \text{bra} \mid \text{ket} \rangle\f$ over \f$\mathbb{R}^D\f$.
 *
 * @tparam D Spatial dimension.
 * @tparam T Scalar type of the bra tree (e.g., `double`, `ComplexDouble`).
 * @tparam U Scalar type of the ket tree.
 * @tparam V Return type deduced as `decltype(T{} * U{})`.
 *
 * @param bra Multiresolution function tree acting as the bra.
 * @param ket Multiresolution function tree acting as the ket.
 * @return The scalar inner product value. For complex types, the bra is
 *         conjugated (i.e., \f$\int \overline{bra(x)}\,ket(x)\,dx\f$).
 *
 * @pre Both trees must be compatible (same MRA/grid conventions).
 * @note Implementations usually reconstruct to consistent representations
 *       before integration; they may refine adaptively to ensure accuracy.
 * @warning For poorly overlapping/aliased grids, the routine may refine
 *          meshes internally, which can be expensive.
 */
template <int D, typename T = double, typename U = T, typename V = decltype(std::declval<T>() * std::declval<U>())>
V dot(FunctionTree<D, T> &bra, FunctionTree<D, U> &ket);

/**
 * @brief Contract two vectors of trees into a scalar field:
 *        \f$out(x) = \sum_i a_i(x)\,\overline{b_i(x)}\f$ (complex) or \f$\sum_i a_i(x)\,b_i(x)\f$ (real).
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient type (`double` or `ComplexDouble`).
 * @param prec Target accuracy for the constructed field (see “Precision semantics” above).
 * @param out  Output scalar field tree receiving the contraction.
 * @param inp_a Vector of factor trees \f$\{a_i\}\f$.
 * @param inp_b Vector of factor trees \f$\{b_i\}\f$; must have the same size and compatible grids as `inp_a`.
 * @param maxIter Maximum refinement passes; `-1` for unlimited.
 * @param absPrec If `true`, interpret `prec` as absolute tolerance.
 *
 * @details
 * Builds the *pointwise* contraction of two equally sized collections of trees,
 * summing products component-wise. This is often used to assemble densities
 * or overlaps distributed over space.
 *
 * The routine adaptively refines `out` to meet `prec`. Input nodes may be
 * transiently reconstructed to compatible representations.
 */
template <int D, typename T>
void dot(double prec, FunctionTree<D, T> &out, FunctionTreeVector<D, T> &inp_a, FunctionTreeVector<D, T> &inp_b, int maxIter = -1, bool absPrec = false);

/**
 * @brief Fast contraction based on node norms (cheap estimate/upper bound).
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient type.
 * @param bra First tree.
 * @param ket Second tree.
 * @param exact If `true`, request exact inner product instead of a norm-based estimate (implementation-dependent).
 * @return A scalar quantity derived from node-wise norms; commonly used
 *         as a quick upper bound or cheap similarity measure.
 *
 * @note When `exact=true`, implementations may fall back to the same
 *       evaluation as @ref dot(bra, ket). If exact evaluation is not
 *       available, `exact` may be ignored.
 */
template <int D, typename T>
double node_norm_dot(FunctionTree<D, T> &bra, FunctionTree<D, T> &ket, bool exact = false);

/**
 * @brief Pointwise product of two trees with a scalar prefactor:
 *        \f$out \leftarrow c \cdot a \cdot (\mathrm{conj}\,b \text{ if requested})\f$.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient type.
 * @param prec Target accuracy for `out`.
 * @param out  Output tree receiving the product.
 * @param c    Scalar prefactor applied to the product.
 * @param inp_a First factor.
 * @param inp_b Second factor.
 * @param maxIter Maximum refinement passes; `-1` for unlimited.
 * @param absPrec If `true`, interpret `prec` as absolute tolerance.
 * @param useMaxNorms If `true`, use max-norm heuristics to guide refinement (may be faster, slightly more conservative).
 * @param conjugate If `true` and `T` is complex, conjugate the **first** factor (bra–ket convention).
 *
 * @details
 * Produces an adaptively refined tree such that the representation error of
 * the pointwise product does not exceed `prec` under the chosen policy.
 */
template <int D, typename T>
void multiply(double prec,
              FunctionTree<D, T> &out,
              T c,
              FunctionTree<D, T> &inp_a,
              FunctionTree<D, T> &inp_b,
              int maxIter = -1,
              bool absPrec = false,
              bool useMaxNorms = false,
              bool conjugate = false);

/**
 * @brief Pointwise product of an arbitrary number of trees:
 *        \f$out \leftarrow \prod_{i} f_i\f$ (optional conjugation of the first factor).
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient type.
 * @param prec Target accuracy for `out`.
 * @param out  Output tree receiving the product.
 * @param inp  List of input tree pointers (non-null, compatible MRAs).
 * @param maxIter Maximum refinement passes; `-1` for unlimited.
 * @param absPrec If `true`, interpret `prec` as absolute tolerance.
 * @param useMaxNorms If `true`, enable max-norm driven refinement.
 * @param conjugate If `true` and `T` is complex, conjugate the **first** factor only.
 *
 * @note The algorithm typically multiplies factors incrementally with
 *       intermediate refinement; ordering can affect performance.
 */
template <int D, typename T>
void multiply(double prec, FunctionTree<D, T> &out, std::vector<FunctionTree<D, T> *> &inp, int maxIter = -1, bool absPrec = false, bool useMaxNorms = false, bool conjugate = false);

/**
 * @brief Pointwise product of a vector of trees:
 *        \f$out \leftarrow \prod_{i} f_i\f$ (optional conjugation of the first factor).
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient type.
 * @param prec Target accuracy for `out`.
 * @param out  Output tree receiving the product.
 * @param inp  Vector wrapper containing input trees (and possibly per-tree scalars).
 * @param maxIter Maximum refinement passes; `-1` for unlimited.
 * @param absPrec If `true`, interpret `prec` as absolute tolerance.
 * @param useMaxNorms If `true`, enable max-norm driven refinement.
 * @param conjugate If `true` and `T` is complex, conjugate the **first** factor only.
 */
template <int D, typename T>
void multiply(double prec, FunctionTree<D, T> &out, FunctionTreeVector<D, T> &inp, int maxIter = -1, bool absPrec = false, bool useMaxNorms = false, bool conjugate = false);

/**
 * @brief Element-wise power: \f$out(x) = \big(inp(x)\big)^{p}\f$.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient type.
 * @param prec Target accuracy for `out`.
 * @param out  Output tree.
 * @param inp  Input tree.
 * @param p    Real exponent.
 * @param maxIter Maximum refinement passes; `-1` for unlimited.
 * @param absPrec If `true`, interpret `prec` as absolute tolerance.
 *
 * @warning For real `T`, negative bases with non-integer `p` are undefined.
 *          For complex `T`, the principal branch is typically used.
 */
template <int D, typename T>
void power(double prec, FunctionTree<D, T> &out, FunctionTree<D, T> &inp, double p, int maxIter = -1, bool absPrec = false);

/**
 * @brief Element-wise square.
 *
 * @tparam D Spatial dimension.
 * @tparam T Coefficient type.
 * @param prec Target accuracy for `out`.
 * @param out  Output tree.
 * @param inp  Input tree.
 * @param maxIter Maximum refinement passes; `-1` for unlimited.
 * @param absPrec If `true`, interpret `prec` as absolute tolerance.
 * @param conjugate If `true` and `T` is complex, compute squared magnitude:
 *                  \f$out = inp \cdot \overline{inp}\f$; otherwise compute
 *                  \f$out = inp \cdot inp\f$.
 */
template <int D, typename T>
void square(double prec, FunctionTree<D, T> &out, FunctionTree<D, T> &inp, int maxIter = -1, bool absPrec = false, bool conjugate = false);

} // namespace mrcpp