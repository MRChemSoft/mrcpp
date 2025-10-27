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

#include <iostream>
#include <vector>

#include "GaussFunc.h"
#include "RepresentableFunction.h"

#include "MRCPP/mrcpp_declarations.h"

namespace mrcpp {

#define GAUSS_EXP_PREC 1.e-10

/** @class GaussExp
 *  @tparam D Spatial dimension (1, 2, 3, …).
 *
 *  @brief Container for a finite linear combination (“expansion”) of
 *  Cartesian Gaussian-type primitives in D dimensions.
 *
 *  Mathematical model
 *  ------------------
 *  - 1D:
 *    \f[
 *      g(x) = \sum_{m=1}^M g_m(x)
 *            = \sum_{m=1}^M \alpha_m \exp\!\big(-\beta_m (x - x_m)^2\big).
 *    \f]
 *  - D dimensions (separable Cartesian form):
 *    \f[
 *      G(\mathbf{x}) = \sum_{m=1}^M G_m(\mathbf{x})
 *                     = \sum_{m=1}^M \prod_{d=1}^D g^{(d)}_m(x_d),
 *    \f]
 *    where each term is represented by a @ref Gaussian<D> (base class) and is
 *    concretely either a pure Gaussian @ref GaussFunc<D> or a Gaussian times a
 *    Cartesian polynomial @ref GaussPoly<D>.
 *
 *  Ownership & invariants
 *  ----------------------
 *  - The expansion OWNS its terms via raw pointers (@c Gaussian<D>*). It
 *    performs deep copies on copy construction/assignment and deletes terms
 *    in the destructor.
 *  - @c funcs[i] is either non-null (a valid Gaussian term) or null for
 *    “empty” slots when constructed with a fixed number of terms.
 *
 *  Typical uses
 *  ------------
 *  - Build analytic functions as sums of Gaussians, evaluate them pointwise
 *    (@ref evalf).
 *  - Combine expansions algebraically: @ref add, @ref mult (by expansion,
 *    single term, polynomial term, or scalar).
 *  - Compute norms and normalize: @ref calcSquareNorm, @ref normalize.
 *  - Manage crude support/visibility via screening: @ref calcScreening,
 *    @ref setScreen.
 */
template <int D> class GaussExp : public RepresentableFunction<D, double> {
public:
    /**
     * @brief Construct an expansion with @p nTerms empty slots.
     *
     * @param nTerms Number of entries reserved in @ref funcs (default 0).
     * @param prec   Historical argument (unused here); kept for API symmetry.
     *
     * After construction, populate terms via @ref setFunc or @ref append.
     */
    GaussExp(int nTerms = 0, double prec = GAUSS_EXP_PREC);

    /** @brief Deep-copy constructor (clones every term via virtual copy()). */
    GaussExp(const GaussExp<D> &gExp);

    /** @brief Deep-copy assignment (existing terms are discarded then cloned). */
    GaussExp &operator=(const GaussExp<D> &gExp);

    /** @brief Destructor: deletes all owned terms and clears pointers. */
    ~GaussExp() override;

    // ---- STL-style iteration over owned pointers (non-const and const) ----
    auto begin() { return funcs.begin(); }
    auto end() { return funcs.end(); }
    const auto begin() const { return funcs.begin(); }
    const auto end() const { return funcs.end(); }

    // ---- Analysis helpers ---------------------------------------------------

    /**
     * @brief Coulomb self-energy of the expansion.
     * @details Implemented for D=3 (see .cpp); throws/not-implemented for others.
     * @note For physical correctness each term should be charge-normalized.
     */
    double calcCoulombEnergy() const;

    /**
     * @brief Compute the squared L2 norm of the expansion:
     *        \f$ \| \sum_i f_i \|_2^2 = \sum_i \|f_i\|^2 + 2\sum_{i<j}\langle f_i,f_j\rangle \f$.
     * @details Uses closed-form Gaussian overlaps (see implementation).
     */
    double calcSquareNorm() const;

    /** @brief Scale all coefficients so that @ref calcSquareNorm() == 1. */
    void normalize();

    /**
     * @brief Precompute per-term screening extents using @p nStdDev “sigma”.
     * @param nStdDev Number of standard deviations used to truncate/screen.
     * @details Also stores the value locally for future reference.
     */
    void calcScreening(double nStdDev = defaultScreening);

    // ---- RepresentableFunction interface ------------------------------------

    /**
     * @brief Pointwise evaluation \f$ f(\mathbf{r}) = \sum_i f_i(\mathbf{r}) \f$.
     * @param r D-dimensional coordinate.
     */
    double evalf(const Coord<D> &r) const override;

    // ---- Other transforms/utilities -----------------------------------------

    /**
     * @brief Build a periodified version of the expansion by tiling each term.
     * @param period Period per axis (e.g., {Lx, Ly, Lz} in 3D).
     * @param nStdDev Screening control for how many images to include.
     */
    GaussExp<D> periodify(const std::array<double, D> &period, double nStdDev = 4.0) const;

    /**
     * @brief Component-wise derivative \f$\partial/\partial x_{\text{dir}}\f$.
     * @param dir Axis index in [0, D).
     * @return New expansion with each term differentiated.
     */
    GaussExp<D> differentiate(int dir) const;

    // ---- Algebra: additive and multiplicative combinators -------------------

    /** @brief Concatenate two expansions (returns all terms from both). */
    GaussExp<D> add(GaussExp<D> &g);

    /** @brief Append a single Gaussian term to this expansion (returns new sum). */
    GaussExp<D> add(Gaussian<D> &g);

    /**
     * @brief Distribute product over terms:
     *        (Σ f_i) * (Σ g_j) = Σ_{ij} f_i⋅g_j (resulting in GaussPoly terms).
     */
    GaussExp<D> mult(GaussExp<D> &g);

    /** @brief Multiply by a single pure Gaussian (resulting in GaussPoly terms). */
    GaussExp<D> mult(GaussFunc<D> &g);

    /** @brief Multiply by a single Gaussian–polynomial term. */
    GaussExp<D> mult(GaussPoly<D> &g);

    /** @brief Return a copy scaled by scalar @p d. */
    GaussExp<D> mult(double d);

    /** @brief Scale coefficients in place by scalar @p d. */
    void multInPlace(double d);

    // ---- Operator sugar (forward to the methods above) ----------------------

    GaussExp<D> operator+(GaussExp<D> &g) { return this->add(g); }
    GaussExp<D> operator+(Gaussian<D> &g) { return this->add(g); }
    GaussExp<D> operator*(GaussExp<D> &g) { return this->mult(g); }
    GaussExp<D> operator*(GaussFunc<D> &g) { return this->mult(g); }
    GaussExp<D> operator*(GaussPoly<D> &g) { return this->mult(g); }
    GaussExp<D> operator*(double d) { return this->mult(d); }
    void operator*=(double d) { this->multInPlace(d); }

    // ---- Accessors ----------------------------------------------------------

    /** @brief Current screening parameter (sign may encode “enabled/disabled”). */
    double getScreening() const { return screening; }

    /** @brief Exponent(s) α per axis for term i. */
    std::array<double, D> getExp(int i) const { return this->funcs[i]->getExp(); }

    /** @brief Scalar coefficient for term i. */
    double getCoef(int i) const { return this->funcs[i]->getCoef(); }

    /** @brief Powers (Cartesian angular momenta) per axis for term i. */
    const std::array<int, D> &getPower(int i) const { return this->funcs[i]->getPower(); }

    /** @brief Center position per axis for term i. */
    const std::array<double, D> &getPos(int i) const { return this->funcs[i]->getPos(); }

    /** @brief Number of (owned) terms in the expansion. */
    int size() const { return this->funcs.size(); }

    /** @brief Mutable access to term i (reference). */
    Gaussian<D> &getFunc(int i) { return *this->funcs[i]; }

    /** @brief Const access to term i (reference). */
    const Gaussian<D> &getFunc(int i) const { return *this->funcs[i]; }

    /** @brief Mutable pointer access (may be null if slot is empty). */
    Gaussian<D> *operator[](int i) { return this->funcs[i]; }

    /** @brief Const pointer access (may be null if slot is empty). */
    const Gaussian<D> *operator[](int i) const { return this->funcs[i]; }

    // ---- Mutators -----------------------------------------------------------

    /**
     * @brief Install a Gaussian–polynomial term at slot i, scaling its coef by c.
     * @details Replaces any existing object at slot i (deletes old).
     */
    void setFunc(int i, const GaussPoly<D> &g, double c = 1.0);

    /**
     * @brief Install a pure Gaussian term at slot i, scaling its coef by c.
     * @details Replaces any existing object at slot i (deletes old).
     */
    void setFunc(int i, const GaussFunc<D> &g, double c = 1.0);

    /**
     * @brief Set global default screening for newly created instances.
     * @throws If @p screen is negative.
     */
    void setDefaultScreening(double screen);

    /**
     * @brief Enable/disable screening for this expansion and forward to terms.
     * @details Conventionally, a positive @ref screening means “enabled” and
     *          a negative value means “disabled”.
     */
    void setScreen(bool screen);

    /** @brief Set (isotropic) exponent(s) α for term i. */
    void setExp(int i, double a) { this->funcs[i]->setExp(a); }

    /** @brief Set scalar coefficient for term i. */
    void setCoef(int i, double b) { this->funcs[i]->setCoef(b); }

    /** @brief Set Cartesian powers for term i. */
    void setPow(int i, const std::array<int, D> &power) { this->funcs[i]->setPow(power); }

    /** @brief Set center position for term i. */
    void setPos(int i, const std::array<double, D> &pos) { this->funcs[i]->setPos(pos); }

    /** @brief Append a single (cloned) Gaussian to the end of the expansion. */
    void append(const Gaussian<D> &g);

    /** @brief Append all terms (cloned) from another expansion. */
    void append(const GaussExp<D> &g);

    /** @brief Stream pretty-printer: prints a summary and the terms. */
    friend std::ostream &operator<<(std::ostream &o, const GaussExp<D> &gExp) { return gExp.print(o); }

    /** @brief Grant @ref Gaussian access to internals where necessary. */
    friend class Gaussian<D>;

protected:
    /** @brief Owned list of Gaussian terms (raw-pointer ownership). */
    std::vector<Gaussian<D> *> funcs;

    /** @brief Default screening parameter for new instances of this @c D. */
    static double defaultScreening;

    /** @brief Instance screening parameter (sign may encode enabled/disabled). */
    double screening{0.0};

    /** @brief Implementation of stream printing (called by operator<<). */
    std::ostream &print(std::ostream &o) const;

    /**
     * @brief Coarse visibility test for adaptive algorithms.
     * @details Returns false if any term declares itself not visible
     *          at the given scale/sample count; true otherwise.
     */
    bool isVisibleAtScale(int scale, int nPts) const override;

    /**
     * @brief Quick zero check on a box \f$[lb,ub]^D\f$:
     *        returns true only if every term is zero on the box.
     */
    bool isZeroOnInterval(const double *lb, const double *ub) const override;
};

} // namespace mrcpp