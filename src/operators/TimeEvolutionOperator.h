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
 * @file TimeEvolutionOperator.h
 * @brief Interface for a separable multiwavelet representation of the
 *        free-particle Schrödinger time-evolution semigroup.
 *
 * The operator approximates (real or imaginary parts of)
 * \f[
 *   U(t) \;=\; e^{\, i\,t\,\Delta}
 * \f]
 * by building an operator tree via cross-correlations between scaling functions
 * and a kernel whose coefficients are expressed through power integrals
 * \f$ \widetilde J_m \f$. Two construction modes are exposed:
 *  - **Uniform** to a user-specified finest scale.
 *  - **Adaptive** down to a fixed scale (bounded work in power integrals).
 *
 * See the .cpp for build details and post-processing steps (MW transform,
 * rough-scale filtering, caching, etc.).
 */

#pragma once

#include "ConvolutionOperator.h"
#include "MWOperator.h"
#include "core/SchrodingerEvolution_CrossCorrelation.h"

namespace mrcpp {

/**
 * @class TimeEvolutionOperator
 * @ingroup operators
 *
 * @brief Multiwavelet operator for the free-particle Schrödinger semigroup.
 *
 * @tparam D Spatial dimensionality (1, 2, or 3).
 *
 * @details
 * Provides a separable @ref ConvolutionOperator-like interface that assembles
 * the matrix elements of
 * \f$ U(t) = e^{\, i\,t\,\Delta} \f$
 * (or its real/imaginary part) in a multi-resolution setting. The actual
 * operator blocks can be accessed via
 * @code
 *   getComponent(0, 0)
 * @endcode
 * after construction (rank-1 expansion in current implementation).
 *
 * Internally, coefficients are generated from per-scale power integrals
 * \f$ \widetilde J_m \f$ and a dedicated cross-correlation calculator suited
 * for the Schrödinger kernel.
 *
 * @note Current implementation targets Legendre scaling functions; practical
 *       use has primarily focused on 1D, but the interface is templated in @p D.
 *
 * @todo Extend to general dimension on arbitrary intervals \f$[a,b]\f$.
 */
template <int D>
class TimeEvolutionOperator : public ConvolutionOperator<D> // One can use ConvolutionOperator instead as well
{
public:
    /**
     * @brief Construct a **uniform** time-evolution operator.
     *
     * @param mra           Target @ref MultiResolutionAnalysis (domain/basis).
     * @param prec          Build precision controlling pruning and tolerances.
     * @param time          Time parameter \f$ t \f$.
     * @param finest_scale  Finest (uniform) scale to which the operator is built.
     * @param imaginary     If `true` build the imaginary part; otherwise real part.
     * @param max_Jpower    Maximum number of power-integral terms (default: 30).
     */
    TimeEvolutionOperator(const MultiResolutionAnalysis<D> &mra,
                          double prec,
                          double time,
                          int finest_scale,
                          bool imaginary,
                          int max_Jpower = 30);

    /**
     * @brief Construct an **adaptive** time-evolution operator.
     *
     * @param mra         Target @ref MultiResolutionAnalysis (domain/basis).
     * @param prec        Build precision controlling pruning and tolerances.
     * @param time        Time parameter \f$ t \f$.
     * @param imaginary   If `true` build the imaginary part; otherwise real part.
     * @param max_Jpower  Maximum number of power-integral terms (default: 30).
     *
     * @details
     * The adaptive build proceeds down to a fixed scale to bound the number of
     * required power integrals; see the source for the current depth choice.
     */
    TimeEvolutionOperator(const MultiResolutionAnalysis<D> &mra,
                          double prec,
                          double time,
                          bool imaginary,
                          int max_Jpower = 30);

    TimeEvolutionOperator(const TimeEvolutionOperator &oper) = delete;            ///< Non-copyable
    TimeEvolutionOperator &operator=(const TimeEvolutionOperator &oper) = delete; ///< Non-assignable
    virtual ~TimeEvolutionOperator() = default;

    /// @return The build precision used to assemble the operator.
    double getBuildPrec() const { return this->build_prec; }

protected:
    /** @name Builder entry points (implementation detail)
     *  Internal construction routines used by the public constructors.
     */
    ///@{
    /**
     * @brief Uniform build to @p finest_scale.
     * @param time         Time parameter \f$ t \f$.
     * @param finest_scale Finest scale to which the tree is constructed.
     * @param imaginary    Build imaginary (true) or real (false) part.
     * @param max_Jpower   Maximum number of power-integral terms.
     */
    void initialize(double time, int finest_scale, bool imaginary, int max_Jpower);

    /**
     * @brief Adaptive build (fixed maximum depth).
     * @param time        Time parameter \f$ t \f$.
     * @param imaginary   Build imaginary (true) or real (false) part.
     * @param max_Jpower  Maximum number of power-integral terms.
     */
    void initialize(double time, bool imaginary, int max_Jpower);

    /**
     * @brief Semi-uniform prototype (not implemented).
     * @param time        Time parameter \f$ t \f$.
     * @param imaginary   Build imaginary (true) or real (false) part.
     * @param max_Jpower  Maximum number of power-integral terms.
     * @warning This method is a placeholder and aborts if called.
     */
    void initializeSemiUniformly(double time, bool imaginary, int max_Jpower);
    ///@}

    /// Set the build precision recorded by this operator.
    void setBuildPrec(double prec) { this->build_prec = prec; }

    double build_prec{-1.0};                                   ///< Build precision (assembly/pruning).
    SchrodingerEvolution_CrossCorrelation *cross_correlation{nullptr}; ///< Per-dimension cross-correlation engine.
};

} // namespace mrcpp