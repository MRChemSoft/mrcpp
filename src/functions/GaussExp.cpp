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
 * @file GaussExp.cpp
 *
 * @brief Implementation of @c GaussExp<D>, a small container for a linear
 *        combination (expansion) of Cartesian Gaussian primitives and/or
 *        Gaussian–polynomial terms. The class offers:
 *        - basic construction/assignment and memory ownership of terms,
 *        - pointwise evaluation,
 *        - algebra (sum/product by distributing over terms),
 *        - norm and normalization helpers,
 *        - crude visibility/screening support,
 *        - Coulomb energy (specialized for D=3),
 *        - periodification helper.
 *
 * Design notes
 * ------------
 * - The expansion holds owning pointers to @c Gaussian<D> (base type), and
 *   concrete terms are either @c GaussFunc<D> (pure Gaussian) or
 *   @c GaussPoly<D> (Gaussian times a Cartesian polynomial).
 * - Operations that combine expansions rely on @c dynamic_cast to handle the
 *   two concrete term types and produce a @c GaussPoly<D> when multiplying.
 * - @b Ownership: this class allocates copies on insert/append and frees them
 *   in the destructor; copy constructor and assignment perform deep copies.
 * - Screening: @c screening is a scalar that configures per-term screening
 *   (e.g., via “n standard deviations”); negative values can be used as a
 *   disabled flag (see @c setScreen). Each term also receives the screen state.
 */

#include "GaussExp.h"

#include <cstdlib>
#include <iostream>

#include "function_utils.h"
#include "GaussFunc.h"
#include "GaussPoly.h"
#include "Gaussian.h"
#include "utils/Printer.h"
#include "utils/math_utils.h"

using namespace Eigen;

namespace mrcpp {

/** @brief Default screening parameter (in “number of standard deviations”).
 *
 * Each dimensional specialization gets its own static. Positive means enabled
 * by default; see @ref setScreen to flip the sign and propagate to terms.
 */
template <int D> double GaussExp<D>::defaultScreening = 10.0;

/**
 * @brief Construct an expansion with a fixed number of (empty) slots.
 *
 * @param nTerms Number of terms (initial capacity).
 * @param prec   Unused here (historical signature compatibility).
 *
 * The vector is filled with @c nullptr placeholders; actual terms must be
 * installed via @ref setFunc or @ref append before use.
 */
template <int D> GaussExp<D>::GaussExp(int nTerms, double /*prec*/) {
    for (int i = 0; i < nTerms; i++) { this->funcs.push_back(nullptr); }
}

/**
 * @brief Deep-copy constructor.
 *
 * Clones each term by calling its virtual @c copy() (polymorphic copy).
 * The @c screening flag/value is copied as well.
 */
template <int D> GaussExp<D>::GaussExp(const GaussExp<D> &gexp) {
    screening = gexp.screening;
    for (unsigned int i = 0; i < gexp.size(); i++) {
        Gaussian<D> *gauss = gexp.funcs[i]->copy();
        this->funcs.push_back(gauss);
    }
}

/**
 * @brief Destructor: deletes all owned terms (if any) and nulls pointers.
 */
template <int D> GaussExp<D>::~GaussExp() {
    for (int i = 0; i < size(); i++) {
        if (this->funcs[i] != nullptr) {
            delete this->funcs[i];
            this->funcs[i] = nullptr;
        }
    }
}

/**
 * @brief Deep-copy assignment (strong exception safety not guaranteed).
 *
 * Existing terms are discarded; the right-hand side is cloned term by term.
 * The @c screening parameter is @b not overwritten (commented line preserves
 * current object’s screening), so only structure/terms are copied.
 */
template <int D> GaussExp<D> &GaussExp<D>::operator=(const GaussExp<D> &gexp) {
    if (&gexp == this) return *this;
    // screening = gexp.screening;
    this->funcs.clear();
    for (unsigned int i = 0; i < gexp.size(); i++) {
        if (gexp.funcs[i] == nullptr) {
            this->funcs.push_back(nullptr);
        } else {
            Gaussian<D> *gauss = gexp.getFunc(i).copy();
            this->funcs.push_back(gauss);
        }
    }
    return *this;
}

/**
 * @brief Pointwise evaluation: sum of all term evaluations at @p r.
 *
 * @param r D-dimensional coordinate.
 * @return  Σ_i term_i(r).
 */
template <int D> double GaussExp<D>::evalf(const Coord<D> &r) const {
    double val = 0.0;
    for (int i = 0; i < this->size(); i++) { val += this->getFunc(i).evalf(r); }
    return val;
}

/**
 * @brief Quick “visibility” test at a given scale and sample count.
 *
 * @details Returns @c false if any term is not visible (fails its own
 * visibility criterion); only if all are visible does it return @c true.
 * This is a conservative conjunction useful for pruning.
 */
template <int D> bool GaussExp<D>::isVisibleAtScale(int scale, int nPts) const {
    for (unsigned int i = 0; i < this->size(); i++) {
        if (not this->getFunc(i).isVisibleAtScale(scale, nPts)) { return false; }
    }
    return true;
}

/**
 * @brief Check whether the expansion is identically zero on [lb,ub]^D.
 *
 * @details Returns @c false if any term says it is non-zero on the box;
 * otherwise returns @c true. Used for quick region elimination.
 */
template <int D> bool GaussExp<D>::isZeroOnInterval(const double *lb, const double *ub) const {
    for (unsigned int i = 0; i < this->size(); i++) {
        if (not this->getFunc(i).isZeroOnInterval(lb, ub)) { return false; }
    }
    return true;
}

/**
 * @brief Install a @c GaussPoly term into slot @p i, scaling its coefficient.
 *
 * @param i Slot index (0-based).
 * @param g Source Gaussian–polynomial term (copied).
 * @param c Extra scalar factor applied multiplicatively to the stored term’s
 *          existing coefficient (so final coef = c * g.coef()).
 */
template <int D> void GaussExp<D>::setFunc(int i, const GaussPoly<D> &g, double c) {
    if (i < 0 or i > (this->size() - 1)) {
        MSG_ERROR("Index out of bounds!");
        return;
    }
    if (this->funcs[i] != nullptr) { delete this->funcs[i]; }
    this->funcs[i] = new GaussPoly<D>(g);
    double coef = this->funcs[i]->getCoef();
    this->funcs[i]->setCoef(c * coef);
}

/**
 * @brief Install a pure @c GaussFunc term into slot @p i, scaling its coefficient.
 *
 * Same semantics as the GaussPoly overload.
 */
template <int D> void GaussExp<D>::setFunc(int i, const GaussFunc<D> &g, double c) {
    if (i < 0 or i > (this->size() - 1)) {
        MSG_ERROR("Index out of bounds!");
        return;
    }
    if (this->funcs[i] != nullptr) { delete this->funcs[i]; }
    this->funcs[i] = new GaussFunc<D>(g);
    double coef = this->funcs[i]->getCoef();
    this->funcs[i]->setCoef(c * coef);
}

/**
 * @brief Append a new term by polymorphic copy.
 */
template <int D> void GaussExp<D>::append(const Gaussian<D> &g) {
    Gaussian<D> *gp = g.copy();
    this->funcs.push_back(gp);
}

/**
 * @brief Append all terms from another expansion (deep copies).
 */
template <int D> void GaussExp<D>::append(const GaussExp<D> &g) {
    for (int i = 0; i < g.size(); i++) {
        Gaussian<D> *gp = g.getFunc(i).copy();
        this->funcs.push_back(gp);
    }
}

/**
 * @brief Differentiate each term with respect to coordinate @p dir and return a new expansion.
 *
 * @param dir Axis index (0..D-1).
 */
template <int D> GaussExp<D> GaussExp<D>::differentiate(int dir) const {
    assert(dir >= 0 and dir < D);
    GaussExp<D> result;
    for (int i = 0; i < this->size(); i++) result.append(this->getFunc(i).differentiate(dir));
    return result;
}

/**
 * @brief Termwise concatenation (sum) with another expansion.
 *
 * @details Produces an expansion containing all terms from @c *this followed
 * by all terms from @p g, by cloning. Coefficients remain unchanged.
 */
template <int D> GaussExp<D> GaussExp<D>::add(GaussExp<D> &g) {
    int nsum = this->size() + g.size();
    GaussExp<D> sum = GaussExp<D>(nsum);

    int n = 0;
    for (int i = 0; i < this->size(); i++) {
        sum.funcs[n] = this->funcs[i]->copy();
        n++;
    }
    for (int i = 0; i < g.size(); i++) {
        sum.funcs[n] = g.funcs[i]->copy();
        n++;
    }

    return sum;
}

/**
 * @brief Concatenate with a single term @p g (at the end).
 */
template <int D> GaussExp<D> GaussExp<D>::add(Gaussian<D> &g) {
    int nsum = this->size() + 1;
    GaussExp<D> sum = GaussExp<D>(nsum);
    for (int n = 0; n < this->size(); n++) { sum.funcs[n] = this->getFunc(n).copy(); }
    sum.funcs[this->size()] = g.copy();
    return sum;
}

/**
 * @brief Product of two expansions by distributivity.
 *
 * @details For each pair of terms, multiply them (Gaussian×Gaussian or
 * Gaussian×GaussPoly) to produce a @c GaussPoly term which is appended to the
 * result. Type dispatch is handled via @c dynamic_cast and throws on unknown
 * runtime types.
 */
template <int D> GaussExp<D> GaussExp<D>::mult(GaussExp<D> &gexp) {
    GaussExp<D> result;
    for (int i = 0; i < this->size(); i++) {
        for (int j = 0; j < gexp.size(); j++) {
            if (auto *f = dynamic_cast<GaussFunc<D> *>(this->funcs[i])) {
                if (auto *g = dynamic_cast<GaussFunc<D> *>(gexp.funcs[j])) {
                    GaussPoly<D> newTerm = (*g) * (*f);
                    result.append(newTerm);
                } else if (auto *g = dynamic_cast<GaussPoly<D> *>(gexp.funcs[j])) {
                    GaussPoly<D> newTerm = (*g) * (*f);
                    result.append(newTerm);
                } else {
                    MSG_ABORT("Invalid Gaussian type!");
                }
            } else if (auto *f = dynamic_cast<GaussPoly<D> *>(this->funcs[i])) {
                if (auto *g = dynamic_cast<GaussFunc<D> *>(gexp.funcs[j])) {
                    GaussPoly<D> newTerm = (*f) * (*g);
                    result.append(newTerm);
                } else if (auto *g = dynamic_cast<GaussPoly<D> *>(gexp.funcs[j])) {
                    GaussPoly<D> newTerm = (*f) * (*g);
                    result.append(newTerm);
                } else {
                    MSG_ABORT("Invalid Gaussian type!");
                }
            } else {
                MSG_ABORT("Invalid Gaussian type!");
            }
        }
    }
    return result;
}

/**
 * @brief Multiply the expansion by a single @c GaussFunc term (distribute over terms).
 */
template <int D> GaussExp<D> GaussExp<D>::mult(GaussFunc<D> &g) {
    GaussExp<D> result;
    int nTerms = this->size();
    for (int n = 0; n < nTerms; n++) {
        if (auto *f = dynamic_cast<GaussFunc<D> *>(this->funcs[n])) {
            GaussPoly<D> newTerm = *f * g;
            result.append(newTerm);
        } else if (auto *f = dynamic_cast<GaussPoly<D> *>(this->funcs[n])) {
            GaussPoly<D> newTerm = *f * g;
            result.append(newTerm);
        } else {
            MSG_ABORT("Invalid Gaussian type!");
        }
    }
    return result;
}

/**
 * @brief Multiply the expansion by a single @c GaussPoly term (distribute over terms).
 */
template <int D> GaussExp<D> GaussExp<D>::mult(GaussPoly<D> &g) {
    int nTerms = this->size();
    GaussExp<D> result(nTerms);
    for (int n = 0; n < nTerms; n++) {
        if (auto *f = dynamic_cast<GaussFunc<D> *>(this->funcs[n])) {
            GaussPoly<D> newTerm(g * *f);
            result.append(newTerm);
        } else if (auto *f = dynamic_cast<GaussPoly<D> *>(this->funcs[n])) {
            GaussPoly<D> newTerm(g * *f);
            result.append(newTerm);
        } else {
            MSG_ABORT("Invalid Gaussian type!");
        }
    }
    return result;
}

/**
 * @brief Return a copy of the expansion scaled by constant @p d.
 */
template <int D> GaussExp<D> GaussExp<D>::mult(double d) {
    GaussExp<D> prod = *this;
    for (int i = 0; i < this->size(); i++) prod.funcs[i]->multConstInPlace(d);
    return prod;
}

/**
 * @brief In-place scaling of all term coefficients by @p d.
 */
template <int D> void GaussExp<D>::multInPlace(double d) {
    for (int i = 0; i < this->size(); i++) this->funcs[i]->multConstInPlace(d);
}

/**
 * @brief Compute \f$\| \sum_i f_i \|_2^2\f$ via self-terms plus cross terms.
 *
 * @details First sum each term’s squared norm, then add the double products
 * (2× overlap) between distinct terms. To ensure closed form overlaps, terms
 * are materialized as @c GaussFunc and @c calcOverlap is used internally.
 */
template <int D> double GaussExp<D>::calcSquareNorm() const {
    /* computing the squares */
    double norm = 0.0;
    for (int i = 0; i < this->size(); i++) {
        double nc = this->funcs[i]->calcSquareNorm();
        norm += nc;
    }
    /* computing the double products */
    for (int i = 0; i < this->size(); i++) {
        GaussExp<D> funcs_i = getFunc(i).asGaussExp(); // Make sure all entries are GaussFunc
        for (int fi = 0; fi < funcs_i.size(); fi++) {
            GaussFunc<D> &func_i = static_cast<GaussFunc<D> &>(funcs_i.getFunc(fi));
            for (int j = i + 1; j < this->size(); j++) {
                GaussExp<D> funcs_j = getFunc(j).asGaussExp(); // Make sure all entries are GaussFunc
                for (int fj = 0; fj < funcs_j.size(); fj++) {
                    GaussFunc<D> &func_j = static_cast<GaussFunc<D> &>(funcs_j.getFunc(fj));
                    double overlap = func_i.calcOverlap(func_j);
                    norm += 2.0 * overlap;
                }
            }
        }
    }
    return norm;
}

/**
 * @brief Normalize the expansion so that @c calcSquareNorm() == 1.
 *
 * @details Scales each term’s coefficient by 1/||f||, where
 * @c ||f|| = sqrt(calcSquareNorm()).
 */
template <int D> void GaussExp<D>::normalize() {
    double norm = std::sqrt(this->calcSquareNorm());
    for (int i = 0; i < this->size(); i++) {
        double coef = this->funcs[i]->getCoef();
        this->funcs[i]->setCoef(coef / norm);
    }
}

/**
 * @brief Set the per-term screening parameter (e.g., n standard deviations).
 *
 * @details Stores @p nStdDev locally and forwards to each term so that they
 * can precompute their own screening envelopes (e.g., bounding radii).
 */
template <int D> void GaussExp<D>::calcScreening(double nStdDev) {
    screening = nStdDev;
    for (int i = 0; i < this->size(); i++) { this->funcs[i]->calcScreening(nStdDev); }
}

/**
 * @brief Enable or disable screening for this expansion and all terms.
 *
 * @param screen If true, make @c screening positive; if false, make it negative.
 *               The sign convention can be used by downstream code as a quick
 *               toggle. Each term receives @c setScreen(screen) as well.
 */
template <int D> void GaussExp<D>::setScreen(bool screen) {
    if (screen) {
        this->screening = std::abs(this->screening);
    } else {
        this->screening = -std::abs(this->screening);
    }
    for (int i = 0; i < this->size(); i++) { this->funcs[i]->setScreen(screen); }
}

// -----------------------------------------------------------------------------
// Project-to-wavelets routine (legacy)
// -----------------------------------------------------------------------------
// The routine below shows how to compute scaling and wavelet coefficients by
// projecting each term separately and expanding to nD via tensor products.
// It is currently commented out (relies on MWNode internals), but the steps
// are left as documentation for future restoration.
/*
template<int D>
void GaussExp<D>::calcWaveletCoefs(MWNode<D> &node) {
    static const int tDim = 1 << D;
    const ScalingBasis &sf = node.getMWTree().getScalingFunctions();
    MatrixXd &scaling = node.getMWTree().getTmpScalingCoefs();
    VectorXd &tmpvec = node.getMWTree().getTmpScalingVector();
    int kp1 = node.getKp1();
    int kp1_d = node.getKp1_d();
    int inpos = kp1_d - kp1;
    int scale = node.getNodeIndex().scale() + 1;
    node.allocCoefs();
    for (int child = 0; child < tDim; child++) {
        int l[D];
        node.calcChildTranslation(child, l);
        for (int n = 0; n < this->size(); n++) {
            if (this->getFunc(n).checkScreen(scale, l)) {
                continue;
            }
            sf.calcScalingCoefs(this->getFunc(n), scale, l, scaling);
            tmpvec.segment(inpos, kp1) = scaling.col(0);
            math_utils::tensorExpandCoefs(D, 0, kp1, kp1_d, scaling, tmpvec);
            node.getCoefs().segment(child * kp1_d, kp1_d) += tmpvec;
        }
    }
    node.mwTransform(Compression);
    node.setHasCoefs();
    node.calcNorms();
}
*/

/**
 * @brief Configure the global default screening parameter for all new instances.
 *
 * @param screen Non-negative value; throws if negative.
 */
template <int D> void GaussExp<D>::setDefaultScreening(double screen) {
    if (screen < 0) { MSG_ERROR("Screening constant cannot be negative!"); }
    defaultScreening = screen;
}

/**
 * @brief Pretty-printer listing the terms (order and parameters).
 */
template <int D> std::ostream &GaussExp<D>::print(std::ostream &o) const {
    o << "Gaussian expansion: " << size() << " terms" << std::endl;
    for (int i = 0; i < size(); i++) {
        o << "Term" << std::setw(3) << i << " :" << std::endl;
        o << getFunc(i) << std::endl << std::endl;
    }
    return o;
}

/**
 * @brief Coulomb self-energy placeholder for general D.
 *
 * @note For D≠3 this is not implemented.
 */
template <int D> double GaussExp<D>::calcCoulombEnergy() const {
    NOT_IMPLEMENTED_ABORT
}

/**
 * @brief Coulomb repulsion energy for D=3 including self-interaction once.
 *
 * @details Loops over pairs (i≤j), expands any composite terms to pure
 * Gaussians, and accumulates @c 2*overlap for i<j and @c 1*overlap for i=j.
 * Each Gaussian is assumed normalized to unit charge
 * \f$ c = (\alpha/\pi)^{3/2} \f$ for physical correctness of the energy.
 */
template <> double GaussExp<3>::calcCoulombEnergy() const {
    double energy = 0.0;
    for (int i = 0; i < this->size(); i++) {
        GaussExp<3> funcs_i = getFunc(i).asGaussExp(); // Make sure all entries are GaussFunc
        for (int fi = 0; fi < funcs_i.size(); fi++) {
            GaussFunc<3> &func_i = static_cast<GaussFunc<3> &>(funcs_i.getFunc(fi));
            for (int j = i; j < this->size(); j++) {
                GaussExp<3> funcs_j = getFunc(j).asGaussExp(); // Make sure all entries are GaussFunc
                for (int fj = 0; fj < funcs_j.size(); fj++) {
                    GaussFunc<3> &func_j = static_cast<GaussFunc<3> &>(funcs_j.getFunc(fj));
                    double c = 2.0;
                    if (i == j) c = 1.0;
                    energy += c * func_i.calcCoulombEnergy(func_j);
                }
            }
        }
    }
    return energy;
}

/**
 * @brief Build a periodified expansion by summing periodic images of each term.
 *
 * @param period  Period vector per axis (Lx, Ly, Lz for D=3).
 * @param nStdDev Controls the width/number of included images (screening).
 * @return A new @c GaussExp whose terms include periodic replicas of the input.
 */
template <int D> GaussExp<D> GaussExp<D>::periodify(const std::array<double, D> &period, double nStdDev) const {
    GaussExp<D> out_exp;
    for (const auto &gauss : *this) {
        auto periodic_gauss = gauss->periodify(period, nStdDev);
        out_exp.append(periodic_gauss);
    }
    return out_exp;
}

// Explicit template instantiations for common dimensions
template class GaussExp<1>;
template class GaussExp<2>;
template class GaussExp<3>;

} // namespace mrcpp