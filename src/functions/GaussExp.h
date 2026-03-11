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

/**
 * @class GaussExp
 * @tparam D Spatial dimension (1, 2, or 3)
 *
 * @brief Gaussian expansion in D dimensions
 *
 * @details
 * - Monodimensional Gaussian expansion:
 *
 * \f$ g(x) = \sum_{m=1}^M g_m(x) = \sum_{m=1}^M \alpha_m e^{-\beta (x-x^0)^2} \f$
 *
 * - Multidimensional Gaussian expansion:
 *
 * \f$ G(x) = \sum_{m=1}^M G_m(x) = \sum_{m=1}^M \prod_{d=1}^D g_m^d(x^d) \f$
 *
 * Each Gaussian-type functions (GTFs) is represented by a @ref Gaussian<D>
 * (base class) and is concretely either a pure Gaussian @ref GaussFunc<D> or a
 * Gaussian times a Cartesian polynomial @ref GaussPoly<D>.
 */
template <int D> class GaussExp : public RepresentableFunction<D, double> {
public:
    /**
     * @brief Construct a Gaussian expansion and initialize each GTF to `nullptr`
     *
     * @param nTerms Number of GTFs (default 0)
     * @param prec Unused here
     *
     * @note After construction, populate GTFs via @ref setFunc or @ref append.
     */
    GaussExp(int nTerms = 0);

    /// @brief Deep-copy constructor (clones every GTF via virtual copy())
    GaussExp(const GaussExp<D> &gExp);

    /// @brief Deep-copy assignment (existing GTFs are discarded then cloned)
    GaussExp &operator=(const GaussExp<D> &gExp);

    /// @brief Destructor: deletes all owned GTFs and clears pointers
    ~GaussExp() override;

    auto begin() { return funcs.begin(); }              ///< @return An iterator pointing to the first GTF
    auto end() { return funcs.end(); }                  ///< @return An iterator pointing to the past-the-end GTF
    const auto begin() const { return funcs.begin(); }  ///< @return A const iterator pointing to the first GTF
    const auto end() const { return funcs.end(); }      ///< @return A const iterator pointing to the past-the-end GTF

    /**
     * @brief Coulomb repulsion energy between all pairs in the Gaussian expansion, including self-interaction
     * @note Each GTF must be normalized to unit charge
     * \f$ c = (\alpha/\pi)^{D/2} \f$ for this to be correct!
     * Currently this function is only implemented for `D=3`.
     */
    double calcCoulombEnergy() const;

    /**
     * @brief Compute the squared L2 norm of the expansion
     * @details Use \f$ \| \sum_i f_i \|_2^2 = \sum_i \|f_i\|^2 + 2\sum_{i<j}\langle f_i,f_j\rangle \f$.
     */
    double calcSquareNorm() const;

    /// @brief Rescale the Gaussian expansion so that its L2 norm equals 1
    void normalize();

    /**
     * @brief Build bounds around the center on each axis and enable screening for each GTF
     * @param nStdDev Number of standard deviations used for the box
     * @note See the implementation of GTF in @ref Gaussian<D>::calcScreening.
     */
    void calcScreening(double nStdDev = defaultScreening);

    /**
     * @brief Evaluate the Gaussian expansion at a D-dimensional coordinate
     * @param r Point (Coord<D>) in physical space in the MRA box
     * @return Gaussian expansion value at the point
     */
    double evalf(const Coord<D> &r) const override;

    /**
     * @brief Generates a Gaussian expansion that is semi-periodic around a unit-cell
     * @param[in] period: The period of the unit cell
     * @param[in] nStdDev: Number of standard diviations covered in each direction. Default 4.0
     * @return Semi-periodic version of a Gaussian expansion around a unit-cell
     *
     * @note See the implementation of each GTF in @ref Gaussian<D>::periodify.
     */
    GaussExp<D> periodify(const std::array<double, D> &period, double nStdDev = 4.0) const;

    /**
     * @brief Analytic derivative d/dx_dir (Cartesian direction) of the Gaussian expansion
     * @param dir Axis index in [0, D-1]
     *
     * @return A GaussExp<D> representing the derivative
     */
    GaussExp<D> differentiate(int dir) const;

    /**
     * @brief Build a new Gaussian expansion that is the combination of this expansion and the other
     * @param g The other Gaussian expansion
     * @return The new Gaussian expansion with all GTFs from this expansion and the other
     */
    GaussExp<D> add(GaussExp<D> &g);

    /**
     * @brief Build a new Gaussian expansion by appending a single GTF to this expansion
     * @param g The single GTF
     * @return The new Gaussian expansion with GTFs from this expansion and the single GTF
     */
    GaussExp<D> add(Gaussian<D> &g);

    /**
     * @brief Build a new Gaussian expansion by multiplying this expansion and the other
     * @param g The other Gaussian expansion
     * @return The new Gaussian expansion \f$ (\sum_{i} f_{i})(\sum_{j} g_{j}) = \sum_{ij} f_{i} g_{j} \f$
     */
    GaussExp<D> mult(GaussExp<D> &g);

    /**
     * @brief Build a new Gaussian expansion by multiplying this expansion and a single GTF
     * @param g The single GTF
     * @return The new Gaussian expansion \f$ (\sum_{i} f_{i}) g = \sum_{i} f_{i} g \f$
     */
    GaussExp<D> mult(GaussFunc<D> &g);

    /**
     * @brief Build a new Gaussian expansion by multiplying this expansion and a single @ref GaussPoly
     * @param g The single @ref GaussPoly
     * @return The new Gaussian expansion \f$ (\sum_{i} f_{i}) g = \sum_{i} g f_{i} \f$
     */
    GaussExp<D> mult(GaussPoly<D> &g);

    /**
     * @brief Build a new Gaussian expansion whose coefficient is scaled by a scalar
     * @param d The scalar
     * @return The new Gaussian expansion */
    GaussExp<D> mult(double d);

    /**
     * @brief Scale coefficients of this expansion in place by a scalar
     * @param d The scalar */
    void multInPlace(double d);

    /**
     * @brief Overload the + operator to return a new Gaussian expansion formed by combining this expansion with the other
     * @param g The other Gaussian expansion
     * @return The new Gaussian expansion \f$ (\sum_{i} f_{i})(\sum_{j} g_{j}) = \sum_{ij} f_{i} g_{j} \f$
     */
    GaussExp<D> operator+(GaussExp<D> &g) { return this->add(g); }
    /**
     * @brief Overload the + operator to return a new Gaussian expansion formed by appending a single GTF to this expansion
     * @param g The single GTF
     * @return The new Gaussian expansion with GTFs from this expansion and the single GTF
     */
    GaussExp<D> operator+(Gaussian<D> &g) { return this->add(g); }
    /**
     * @brief Overload the * operator to return a new Gaussian expansion formed by multiplying this expansion and the other
     * @param g The other Gaussian expansion
     * @return The new Gaussian expansion \f$ (\sum_{i} f_{i})(\sum_{j} g_{j}) = \sum_{ij} f_{i} g_{j} \f$
     */
    GaussExp<D> operator*(GaussExp<D> &g) { return this->mult(g); }
    /**
     * @brief Overload the * operator to return a new Gaussian expansion formed by multiplying this expansion and a single GTF
     * @param g The single GTF
     * @return The new Gaussian expansion \f$ (\sum_{i} f_{i}) g = \sum_{i} f_{i} g \f$
     */
    GaussExp<D> operator*(GaussFunc<D> &g) { return this->mult(g); }
    /**
     * @brief Overload the * operator to return a new Gaussian expansion formed by multiplying this expansion and a single @ref GaussPoly
     * @param g The single @ref GaussPoly
     * @return The new Gaussian expansion \f$ (\sum_{i} f_{i}) g = \sum_{i} g f_{i} \f$
     */
    GaussExp<D> operator*(GaussPoly<D> &g) { return this->mult(g); }
    /**
     * @brief Overload the * operator to return a new Gaussian expansion whose coefficient is scaled by a scalar
     * @param d The scalar
     * @return The new Gaussian expansion */
    GaussExp<D> operator*(double d) { return this->mult(d); }
    /**
     * @brief Overload the * operator to scale coefficients of this expansion in place by a scalar
     * @param d The scalar */
    void operator*=(double d) { this->multInPlace(d); }

    /// @brief Get screening parameter
    double getScreening() const { return screening; }

    /// @brief Get monomial exponent on the axis for the i-th GTF
    std::array<double, D> getExp(int i) const { return this->funcs[i]->getExp(); }

    /// @brief Get coefficient for the i-th GTF
    double getCoef(int i) const { return this->funcs[i]->getCoef(); }

    /// @brief Get powers for the i-th GTF
    const std::array<int, D> &getPower(int i) const { return this->funcs[i]->getPower(); }

    /// @brief Get position for the i-th GTF
    const std::array<double, D> &getPos(int i) const { return this->funcs[i]->getPos(); }

    /// @brief Get number of GTFs in the expansion
    int size() const { return this->funcs.size(); }

    /// @brief Get mutable access to the i-th GTF
    Gaussian<D> &getFunc(int i) { return *this->funcs[i]; }

    /// @brief Get const access to the i-th GTF
    const Gaussian<D> &getFunc(int i) const { return *this->funcs[i]; }

    /// @brief Get mutable pointer access to the i-th GTF
    Gaussian<D> *operator[](int i) { return this->funcs[i]; }

    /// @brief Get const pointer access to the i-th GTF
    const Gaussian<D> *operator[](int i) const { return this->funcs[i]; }

    /**
     * @brief Set a @ref GaussPoly for the i-th GTF in the expansion and scale its coefficient by a scalar
     * @param i The i-th GTF
     * @param g The @ref GaussPoly
     * @param c The scalar
     * @note Existing i-th GTF will be deleted
     */
    void setFunc(int i, const GaussPoly<D> &g, double c = 1.0);

    /**
     * @brief Set a single GTF for the i-th GTF in the expansion and scale its coefficient by a scalar
     * @param i The i-th GTF
     * @param g The @ref GaussPoly
     * @param c The scalar
     * @note Existing i-th GTF will be deleted
     */
    void setFunc(int i, const GaussFunc<D> &g, double c = 1.0);

    /// @brief Set global default screening for the Gaussian expansion
    void setDefaultScreening(double screen);

    /**
     * @brief Enable/disable screening for this expansion and forward to all GTFs
     * @details Conventionally, a positive @ref screening means "enabled" and
     * a negative value means "disabled".
     */
    void setScreen(bool screen);

    /**
     * @brief Set (isotropic) exponent the i-th GTF
     * @param i The i-th GTF
     * @param a The (isotropic) exponent
     */
    void setExp(int i, double a) { this->funcs[i]->setExp(a); }

    /**
     * @brief Set coefficient for the i-th GTF
     * @param i The i-th GTF
     * @param b The coefficient
     */
    void setCoef(int i, double b) { this->funcs[i]->setCoef(b); }

    /**
     * @brief Set powers for the i-th GTF
     * @param i The i-th GTF
     * @param power The powers
     */
    void setPow(int i, const std::array<int, D> &power) { this->funcs[i]->setPow(power); }

    /**
     * @brief Set center coordinates for the i-th GTF
     * @param i The i-th GTF
     * @param pos The center coordinates 
     */
    void setPos(int i, const std::array<double, D> &pos) { this->funcs[i]->setPos(pos); }

    /**
     * @brief Append a single GTF to the end of the expansion
     * @param g The single GTF
     */
    void append(const Gaussian<D> &g);

    /**
     * @brief Append all GTFs from the other expansion
     * @param g The other expansion
     */
    void append(const GaussExp<D> &g);

    /** @brief Stream pretty-printer (delegates to protected function @ref GaussExp<D>::print) */
    friend std::ostream &operator<<(std::ostream &o, const GaussExp<D> &gExp) { return gExp.print(o); }

    friend class Gaussian<D>;

protected:
    std::vector<Gaussian<D> *> funcs;

    static double defaultScreening;

    double screening{0.0};

    /**
     * @brief Implementation of stream printing (called by operator<<)
     * @param o The output stream
     * @return The output stream
     */
    std::ostream &print(std::ostream &o) const;

    /**
     * @brief Heuristic visibility vs. resolution scale and quadrature sampling
     * @param scale Dyadic scale (tile size ~ 2^{-scale})
     * @param nPts Number of quadrature points per tile edge
     * @return false if any GTF declares itself not visible, true otherwise
     */
    bool isVisibleAtScale(int scale, int nPts) const override;

    /**
     * @brief Quick check whether the expansion is essentially zero on [la,lb] per axis
     * @param la Lower bounds array of length D
     * @param øb Upper bounds array of length D
     * @return true only if each GTF is effectively zero on [la,lb]
     */
    bool isZeroOnInterval(const double *lb, const double *ub) const override;
};

} // namespace mrcpp
