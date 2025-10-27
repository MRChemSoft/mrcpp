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
 *  Implementation notes (high-level)
 *  ---------------------------------
 *  This file implements generic (templated) functionality shared by all
 *  Gaussian-like primitives in D dimensions, i.e. @ref Gaussian<D>.
 *
 *  Key responsibilities:
 *    - Store and initialize the parameters of a Cartesian Gaussian:
 *        * alpha[d]  : per-axis exponents  β_d  (>0)
 *        * coef      : scalar prefactor    α
 *        * pos[d]    : center coordinates  R_d
 *        * power[d]  : polynomial powers   p_d ∈ {0,1,2,...}
 *    - Compose two Gaussians into a *pure* Gaussian by completing the square
 *      (multPureGauss), leaving the polynomial factors to higher layers.
 *    - Build cheap screening boxes / visibility tests to avoid unnecessary
 *      work when projecting to grids/trees (calcScreening, checkScreen,
 *      isVisibleAtScale, isZeroOnInterval).
 *    - Provide utility evaluations on batches of points (evalf over matrices).
 *    - Compute overlaps by expanding (if needed) into GaussFunc terms and
 *      using the Obara–Saika 1D recurrences (via function_utils).
 *    - Create semi-periodic images of a Gaussian inside a unit cell (periodify).
 */

#include <numeric>

#include "Gaussian.h"
#include "GaussExp.h"
#include "function_utils.h"
#include "trees/NodeIndex.h"
#include "utils/Printer.h"
#include "utils/details.h"
#include "utils/math_utils.h"

using namespace Eigen;

namespace mrcpp {

/*---------------------------*
 * Constructors / init state *
 *---------------------------*/

/** @brief Isotropic-constructor: fill all D exponents with the same value @p a. */
template <int D>
Gaussian<D>::Gaussian(double a, double c, const Coord<D> &r, const std::array<int, D> &p)
        : screen(false)          // screening disabled by default
        , coef(c)                // scalar amplitude
        , power(p)               // Cartesian powers
        , pos(r) {               // center
    this->alpha.fill(a);         // isotropic exponent β_d = a ∀ d
}

/** @brief Anisotropic-constructor: exponents are provided per axis. */
template <int D>
Gaussian<D>::Gaussian(const std::array<double, D> &a, double c, const Coord<D> &r, const std::array<int, D> &p)
        : screen(false)
        , coef(c)
        , power(p)
        , alpha(a)
        , pos(r) {}

/*----------------------------------------------------*
 * Multiply two *pure* Gaussians (no polynomial part) *
 *----------------------------------------------------*/
/**
 * @brief Complete-the-square product of two Gaussians into this object.
 *
 * Given
 *   G_L(x) = exp[-Σ_d α_L(d) (x_d - R_L(d))^2],
 *   G_R(x) = exp[-Σ_d α_R(d) (x_d - R_R(d))^2],
 * their product is a *single* Gaussian
 *   G_P(x) = C · exp[-Σ_d α_P(d) (x_d - R_P(d))^2],
 * where
 *   α_P(d) = α_L(d) + α_R(d),
 *   μ(d)   = α_L(d) α_R(d) / α_P(d)         (reduced exponent),
 *   R_P(d) = [α_L(d) R_L(d) + α_R(d) R_R(d)] / α_P(d),
 *   C      = exp[-Σ_d μ(d) (R_L(d) - R_R(d))^2].
 *
 * The polynomial prefactors (if any) are handled elsewhere (e.g. GaussFunc→GaussPoly).
 */
template <int D> void Gaussian<D>::multPureGauss(const Gaussian<D> &lhs, const Gaussian<D> &rhs) {

    auto newAlpha = std::array<double, D>{};
    auto mju = std::array<double, D>{};
    for (auto d = 0; d < D; d++) {
        newAlpha[d] = lhs.alpha[d] + rhs.alpha[d];            // α_P = α_L + α_R
        mju[d] = (lhs.alpha[d] * rhs.alpha[d]) / newAlpha[d]; // μ = α_L α_R / (α_L + α_R)
    }
    auto newPos = std::array<double, D>{};
    auto relPos = std::array<double, D>{};

    double newCoef = 1.0;
    for (int d = 0; d < D; d++) {
        // Center of the product (weighted by exponents)
        newPos[d] = (lhs.alpha[d] * lhs.pos[d] + rhs.alpha[d] * rhs.pos[d]) / newAlpha[d];
        relPos[d] = lhs.pos[d] - rhs.pos[d]; // R_L - R_R
        // Normalization factor from completing the square
        newCoef *= std::exp(-mju[d] * std::pow(relPos[d], 2.0));
    }
    setExp(newAlpha);
    setPos(newPos);
    setCoef(newCoef);
}

/*--------------------------------------------*
 * Screening boxes and quick-visibility tests *
 *--------------------------------------------*/
/**
 * @brief Build an axis-aligned bounding box [A,B] that captures
 *        ~erf coverage based on nStdDev standard deviations.
 *
 * For each dimension d, the 1D Gaussian has variance σ_d^2 = 1/(2 α_d).
 * We choose bounds R_d ± nStdDev * σ_d. Setting @c screen=true enables
 * fast culling in eval and tree projection.
 */
template <int D> void Gaussian<D>::calcScreening(double nStdDev) {
    assert(nStdDev > 0);
    if (not this->isBounded()) {
        // Lazy-allocate bounds arrays if needed
        this->bounded = true;
        this->A = new double[D];
        this->B = new double[D];
    }
    for (int d = 0; d < D; d++) {
        double limit = std::sqrt(nStdDev / this->alpha[d]); // nStdDev * σ_d where σ_d = 1/sqrt(2 α_d)
        this->A[d] = this->pos[d] - limit;
        this->B[d] = this->pos[d] + limit;
    }
    screen = true;
}

/**
 * @brief Tile-level screening: test whether a node box [a,b] at scale n
 *        lies entirely outside this Gaussian’s screening box.
 *
 * The physical length of a dyadic tile at scale n is 2^{-n}. The tile's
 * coordinate bounds are computed from its integer translations l[d].
 * If the tile is completely outside [A,B] on any axis, return true
 * (i.e., we can skip processing that tile).
 */
template <int D> bool Gaussian<D>::checkScreen(int n, const int *l) const {
    if (not getScreen()) { return false; }
    double length = std::pow(2.0, -n); // tile size
    const double *A = this->getLowerBounds();
    const double *B = this->getUpperBounds();
    for (int d = 0; d < D; d++) {
        double a = length * l[d];         // tile lower bound in dim d
        double b = length * (l[d] + 1);   // tile upper bound in dim d
        if (a > B[d] or b < A[d]) { return true; } // entirely outside -> culled
    }
    return false;
}

/**
 * @brief Heuristic visibility test vs. resolution scale and quadrature count.
 *
 * A Gaussian of standard deviation σ should not be represented at
 * resolutions finer than ~σ (no additional information). We compare the
 * scale against a heuristic derived from σ and the number of quadrature points.
 */
template <int D> bool Gaussian<D>::isVisibleAtScale(int scale, int nQuadPts) const {
    for (auto &alp : this->alpha) {
        double stdDeviation = std::pow(2.0 * alp, -0.5);                      // σ = 1/√(2α)
        auto visibleScale = static_cast<int>(-std::floor(std::log2(nQuadPts * 0.5 * stdDeviation)));
        // If requested 'scale' is *finer* (smaller length) than what this σ supports, hide it.
        if (scale < visibleScale) return false;
    }

    return true;
}

/**
 * @brief Quick zero test on an interval: returns true if the Gaussian
 *        is negligible on [a,b] (component-wise), using ±5σ rule.
 *
 * If the interval is completely outside [R-5σ, R+5σ] on any axis, we
 * consider the function zero there for practical purposes.
 */
template <int D> bool Gaussian<D>::isZeroOnInterval(const double *a, const double *b) const {
    for (int i = 0; i < D; i++) {
        double stdDeviation = std::pow(2.0 * this->alpha[i], -0.5);
        double gaussBoxMin = this->pos[i] - 5.0 * stdDeviation;
        double gaussBoxMax = this->pos[i] + 5.0 * stdDeviation;
        if (a[i] > gaussBoxMax or b[i] < gaussBoxMin) return true;
    }
    return false;
}

/*---------------------------------------------*
 * Batch evaluation (matrix of points → values) *
 *---------------------------------------------*/
/**
 * @brief Evaluate the *separable* 1D factors on a batch of points.
 *
 * @param[in]  points  Matrix (N×D) of coordinates; column d contains the d-th coordinate of all N points.
 * @param[out] values  Matrix (N×D) to be filled with per-axis factors:
 *                     values(i,d) = g_d( points(i,d) ).
 *
 * Note: this does not multiply across dimensions; higher-level code can
 *       combine the columns (e.g., by product) if the full D-D value is needed.
 */
template <int D> void Gaussian<D>::evalf(const MatrixXd &points, MatrixXd &values) const {
    assert(points.cols() == D);
    assert(points.cols() == values.cols());
    assert(points.rows() == values.rows());
    for (int d = 0; d < D; d++) {
        for (int i = 0; i < points.rows(); i++) { values(i, d) = evalf1D(points(i, d), d); }
    }
}

/*--------------------------------------*
 * Convenience: maximum standard dev σ  *
 *--------------------------------------*/
/**
 * @brief Return the maximum standard deviation across axes.
 *
 * For isotropic exponents, that is 1/√(2 α). For anisotropic, compute
 * σ_d = 1/√(2 α_d) per axis and return max_d σ_d. Used in periodification.
 */
template <int D> double Gaussian<D>::getMaximumStandardDiviation() const {

    if (details::are_all_equal<D>(this->getExp())) {
        auto exponent = this->getExp()[0];
        return 1.0 / std::sqrt(2.0 * exponent);
    } else {
        auto exponents = this->getExp();
        auto standard_deviations = std::array<double, D>{};
        for (auto i = 0; i < D; i++) { standard_deviations[i] = 1.0 / std::sqrt(2.0 * exponents[i]); }
        return *std::max_element(standard_deviations.begin(), standard_deviations.end());
    }
}

/*-------------------------*
 * Overlap ⟨G|G'⟩ utilities *
 *-------------------------*/
/**
 * @brief General overlap by expanding both sides into @ref GaussFunc terms
 *        (if needed) and summing pairwise 1D Obara–Saika products.
 *
 * The helper function_utils::calc_overlap(GaussFunc,GaussFunc) performs the
 * per-dimension recursion and multiplies contributions across D.
 */
template <int D> double Gaussian<D>::calcOverlap(const Gaussian<D> &inp) const {
    const auto &bra_exp = this->asGaussExp(); // Make sure all entries are GaussFunc
    const auto &ket_exp = inp.asGaussExp();   // Make sure all entries are GaussFunc

    double S = 0.0;
    for (int i = 0; i < bra_exp.size(); i++) {
        const auto &bra_i = static_cast<const GaussFunc<D> &>(bra_exp.getFunc(i));
        for (int j = 0; j < ket_exp.size(); j++) {
            const auto &ket_j = static_cast<const GaussFunc<D> &>(ket_exp.getFunc(j));
            S += function_utils::calc_overlap(bra_i, ket_j);
        }
    }
    return S;
}

/*-----------------------------*
 * Semi-periodic “image” clone *
 *-----------------------------*/
/**
 * @brief Build a semi-periodic expansion by replicating this Gaussian on a
 *        Cartesian lattice so that most of the mass (≈erf coverage) lies
 *        within a single unit cell.
 *
 * @param period   Period vector (cell lengths per axis).
 * @param nStdDev  Number of σ to keep around the central copy (default 4.0).
 * @returns        A @ref GaussExp<D> consisting of translated copies.
 *
 * Algorithm:
 *   1) Fold the original center into the principal cell [0,period).
 *   2) Estimate the number of neighbor cells needed so that ±nStdDev·σ fits.
 *   3) Generate all translations in the (2N+1)^D cube around the folded center.
 *   4) Copy and shift the Gaussian for each translation and append to result.
 */
template <int D> GaussExp<D> Gaussian<D>::periodify(const std::array<double, D> &period, double nStdDev) const {
    GaussExp<D> gauss_exp;
    auto pos_vec = std::vector<Coord<D>>();

    auto x_std = nStdDev * this->getMaximumStandardDiviation();

    // This lambda computes how many neighbor cells are needed (per axis)
    // so that the ±x_std window is covered by translated images.
    auto neighbooring_cells = [period, x_std](auto pos) {
        auto needed_cells_vec = std::vector<int>();
        for (auto i = 0; i < D; i++) {
            auto upper_bound = pos[i] + x_std;
            auto lower_bound = pos[i] - x_std;
            (void)lower_bound; // not used explicitly; retained for clarity
            // Minimal number of positive cell steps so that [pos-x_std, pos+x_std] is inside coverage.
            needed_cells_vec.push_back(std::ceil(upper_bound / period[i]));
        }

        return *std::max_element(needed_cells_vec.begin(), needed_cells_vec.end());
    };

    // Fold starting position into the principal cell
    auto startpos = this->getPos();

    for (auto d = 0; d < D; d++) {
        startpos[d] = std::fmod(startpos[d], period[d]);
        if (startpos[d] < 0) startpos[d] += period[d];
    }

    // Symmetric image range: from -N to +N cells in each dimension
    auto nr_cells_upp_and_down = neighbooring_cells(startpos);
    for (auto d = 0; d < D; d++) { startpos[d] -= nr_cells_upp_and_down * period[d]; }

    // Generate a (2N+1)^D Cartesian product of offsets
    auto tmp_pos = startpos;
    std::vector<double> v(2 * nr_cells_upp_and_down + 1);
    std::iota(v.begin(), v.end(), 0.0);
    auto cart = math_utils::cartesian_product(v, D);
    for (auto &c : cart) {
        for (auto i = 0; i < D; i++) c[i] *= period[i];
    }
    // Shift coordinates by the starting corner
    for (auto &c : cart) std::transform(c.begin(), c.end(), tmp_pos.begin(), c.begin(), std::plus<double>());
    // Convert vectors to mrcpp::Coord
    for (auto &c : cart) {
        mrcpp::Coord<D> pos;
        std::copy_n(c.begin(), D, pos.begin());
        pos_vec.push_back(pos);
    }

    // Create the translated copies
    for (auto &pos : pos_vec) {
        auto *gauss = this->copy();
        gauss->setPos(pos);
        gauss_exp.append(*gauss);
        delete gauss;
    }

    return gauss_exp;
}

/*-----------------------------*
 * Explicit template instances *
 *-----------------------------*/
template class Gaussian<1>;
template class Gaussian<2>;
template class Gaussian<3>;

} // namespace mrcpp
