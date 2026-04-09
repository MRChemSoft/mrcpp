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
 * @brief Plotting utilities for MRCPP functions and trees.
 *
 * This header declares a lightweight plotting helper that samples
 * multivariate functions (or visualizes trees) on simple, equidistant
 * grids derived from user-provided span vectors. It supports 1D (line),
 * 2D (surface), and 3D (cube) outputs and can also dump tree grids.
 *
 * The plotting domain is parameterized by an origin @p O and up to three
 * span vectors @p A, @p B, @p C (not required to be orthogonal). For an
 * overview of the sampling conventions, see @ref mrcpp::Plotter.
 */

#include <Eigen/Core>

#include <fstream>
#include <iostream>
#include <map>
#include <string>

#include "MRCPP/mrcpp_declarations.h"

namespace mrcpp {

/**
 * @class Plotter
 * @tparam D Spatial dimension of the *function* being sampled (1–3).
 * @tparam T Scalar type of the function values (e.g., double, ComplexDouble).
 *
 * @brief Sample multivariate functions on equidistant grids and write results.
 *
 * ### Domain definition
 * The sampling region is specified by:
 * - Origin **O**
 * - Span vectors **A**, **B**, **C**
 *
 * The actual plot type determines how these are used:
 * - **Line plot**: points along **A** starting at **O**
 * - **Surface plot**: a 2D lattice in the parallelogram spanned by **A** and **B**
 * - **Cube plot**: a 3D lattice in the parallelotope spanned by **A**, **B**, **C**
 *
 * None of **A**, **B**, **C** need to be orthogonal.
 *
 * ### Output
 * This class writes simple text files (one value per line or a cube-like block)
 * suitable for quick inspection or feeding into downstream visualization tools.
 * Grid export for trees writes a mesh for node boxes (D=3).
 *
 * @note The template parameter @p D reflects the *intrinsic* dimensionality of
 * the function/tree. A 3D function can still be sampled along a 1D line using
 * @ref linePlot by providing only **A** (and leaving **B**, **C** unused).
 */
template <int D, typename T> class Plotter {
public:
    /**
     * @brief Construct a plotter with a given origin.
     * @param o Plot origin (defaults to the zero vector).
     */
    explicit Plotter(const Coord<D> &o = {});
    virtual ~Plotter() = default;

    /**
     * @brief Set the filename suffix for a plot type.
     *
     * @param t Plot type key (see @ref type).
     * @param s Suffix including the dot (e.g., ".line", ".surf").
     *
     * @details The suffix is appended to the base filename passed to the
     * plotting routines. Defaults are set in the constructor.
     */
    void setSuffix(int t, const std::string &s);

    /**
     * @brief Set the plot origin.
     * @param o New origin **O**.
     */
    void setOrigin(const Coord<D> &o);

    /**
     * @brief Define (or update) the plot span vectors.
     *
     * @param a Vector **A** (required).
     * @param b Vector **B** (optional; used for 2D/3D sampling).
     * @param c Vector **C** (optional; used for 3D sampling).
     *
     * @note Vectors are not required to be orthogonal. The number of points
     * per span is given at call time for each plot type.
     */
    void setRange(const Coord<D> &a, const Coord<D> &b = {}, const Coord<D> &c = {});

    /**
     * @brief Write a grid visualization of a function tree.
     *
     * @param tree Multiresolution tree to visualize.
     * @param fname Base filename (suffix for @ref Grid is appended).
     *
     * @details Exports the end-node grid (and roots) of @p tree. The concrete
     * output is implementation-dependent; for D=3 it is a geomview-friendly
     * mesh (see the .grid writer in the implementation).
     *
     * @warning Meaningful only when the implementation supports the given @p D.
     */
    void gridPlot(const MWTree<D, T> &tree, const std::string &fname);

    /**
     * @brief Sample a function along a line @f$ O + s\,A @f$.
     *
     * @param npts Number of equidistant points along **A**; use @c {N}.
     * @param func Function to evaluate.
     * @param fname Base filename (suffix for @ref Line is appended).
     *
     * @details Generates @c npts[0] positions:
     * @f$ r_i = O + \frac{i}{N-1} A,\ i=0,\dots,N-1 @f$
     * and writes coordinates and values in text form.
     *
     * @pre @ref setRange must have set a non-zero **A**; otherwise this call
     * will fail validation.
     */
    void linePlot(const std::array<int, 1> &npts, const RepresentableFunction<D, T> &func, const std::string &fname);

    /**
     * @brief Sample a function on a surface spanned by **A**, **B**.
     *
     * @param npts Number of points along {**A**, **B**}; use @c {Na, Nb}.
     * @param func Function to evaluate.
     * @param fname Base filename (suffix for @ref Surface is appended).
     *
     * @details Generates positions
     * @f$ r_{ij} = O + \frac{i}{N_a-1}A + \frac{j}{N_b-1}B @f$ and writes
     * coordinates and values in text form.
     *
     * @pre @ref setRange must have set non-zero **A** and **B** when used in 2D/3D.
     */
    void surfPlot(const std::array<int, 2> &npts, const RepresentableFunction<D, T> &func, const std::string &fname);

    /**
     * @brief Sample a function in a 3D block spanned by **A**, **B**, **C**.
     *
     * @param npts Number of points along {**A**, **B**, **C**}; use @c {Na, Nb, Nc}.
     * @param func Function to evaluate.
     * @param fname Base filename (suffix for @ref Cube is appended).
     *
     * @details Generates positions
     * @f$ r_{ijk} = O + \frac{i}{N_a-1}A + \frac{j}{N_b-1}B + \frac{k}{N_c-1}C @f$
     * and writes values in a simple cube-like format suitable for volumetric viewers.
     *
     * @pre @ref setRange must have set non-zero **A**, **B**, **C** when used in 3D.
     */
    void cubePlot(const std::array<int, 3> &npts, const RepresentableFunction<D, T> &func, const std::string &fname);

    /**
     * @brief Plot type selector used for file suffix mapping.
     */
    enum type { Line,    /**< 1D sampling along **A** */
                Surface, /**< 2D sampling on **A**–**B** lattice */
                Cube,    /**< 3D sampling on **A**–**B**–**C** lattice */
                Grid     /**< Grid/mesh export for trees */
    };

protected:
    /** @name Plot domain and output state */
    ///@{
    Coord<D> O{}; ///< Plot origin.
    Coord<D> A{}; ///< Span vector for line plots and first lattice axis.
    Coord<D> B{}; ///< Span vector for surface/cube plots (second axis).
    Coord<D> C{}; ///< Span vector for cube plots (third axis).
    std::ofstream fstrm{};           ///< Owned output stream storage.
    std::ofstream *fout{nullptr};    ///< Active output stream (points to @ref fstrm).
    std::map<int, std::string> suffix{}; ///< Per-type filename suffix map.
    ///@}

    /**
     * @brief Compute step size to place @p pts samples along a span.
     * @param vec Span vector (**A**, **B**, or **C**).
     * @param pts Number of points along the span (>= 1).
     * @return Component-wise step equals @f$ \frac{\text{vec}}{\max(1, pts-1)} @f$.
     *
     * @note When @p pts == 1 the single sample is placed at the origin offset,
     * and the step is unused (implementation guards against division by zero).
     */
    Coord<D> calcStep(const Coord<D> &vec, int pts) const;

    /**
     * @brief Generate equidistant coordinates for a line plot.
     * @param pts_a Points along **A**.
     * @return Matrix of size (pts_a × D) with row-wise coordinates.
     */
    Eigen::MatrixXd calcLineCoordinates(int pts_a) const;

    /**
     * @brief Generate equidistant coordinates for a surface plot.
     * @param pts_a Points along **A**.
     * @param pts_b Points along **B**.
     * @return Matrix of size ((pts_a*pts_b) × D) with row-wise coordinates.
     */
    Eigen::MatrixXd calcSurfCoordinates(int pts_a, int pts_b) const;

    /**
     * @brief Generate equidistant coordinates for a cube plot.
     * @param pts_a Points along **A**.
     * @param pts_b Points along **B**.
     * @param pts_c Points along **C**.
     * @return Matrix of size ((pts_a*pts_b*pts_c) × D) with row-wise coordinates.
     */
    Eigen::MatrixXd calcCubeCoordinates(int pts_a, int pts_b, int pts_c) const;

    /**
     * @brief Evaluate a representable function on given coordinates.
     * @param func Function to sample.
     * @param coords Row-major matrix of coordinates (N × D).
     * @return Column vector of values (size N).
     *
     * @note The implementation may use parallel evaluation (e.g., OpenMP)
     * when available at build time.
     */
    Eigen::Matrix<T, Eigen::Dynamic, 1> evaluateFunction(const RepresentableFunction<D, T> &func, const Eigen::MatrixXd &coords) const;

    /**
     * @brief Write coordinates and values as text rows.
     * @param coords Row-major coordinates (N × D).
     * @param values Column vector (size N).
     *
     * @details Each output line contains D coordinates followed by the
     * function value. Floating-point formatting is implementation-defined.
     */
    void writeData(const Eigen::MatrixXd &coords, const Eigen::Matrix<T, Eigen::Dynamic, 1> &values);

    /**
     * @brief Write volumetric (cube) data.
     * @param npts Lattice sizes {Na, Nb, Nc}.
     * @param values Column vector of length Na*Nb*Nc.
     *
     * @details Default implementation may be a stub; specialized versions
     * (e.g., D=3) provide actual volume exporters (cube/voxel formats).
     */
    virtual void writeCube(const std::array<int, 3> &npts, const Eigen::Matrix<T, Eigen::Dynamic, 1> &values);

    /**
     * @brief Write a grid/mesh representation of a tree.
     * @param tree Tree whose node boxes should be exported.
     *
     * @details Implementation targets D=3 (geomview-like mesh). Other
     * dimensionalities may provide no-ops or alternative encodings.
     */
    void writeGrid(const MWTree<D, T> &tree);

    /**
     * @brief Emit a single node's box edges/faces to the active stream.
     * @param node Node to visualize.
     * @param color Renderer-dependent color string (implementation-defined).
     */
    void writeNodeGrid(const MWNode<D, T> &node, const std::string &color);

private:
    /**
     * @brief Validate that required span vectors are non-zero.
     * @param dim Required plot dimensionality (1, 2, or 3).
     * @return @c true if all needed spans (**A**, **B**, **C**) are non-zero.
     */
    bool verifyRange(int dim) const;

    /**
     * @brief Open/prepare the output file stream.
     * @param fname Base filename plus suffix (if non-empty).
     *
     * @details If @p fname is empty, reuses the current stream; otherwise
     * closes any previous stream and opens the new one.
     */
    void openPlot(const std::string &fname);

    /**
     * @brief Close the output stream if open and reset state.
     */
    void closePlot();
};

} // namespace mrcpp