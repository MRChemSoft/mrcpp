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

#include <Eigen/Core>

#include <fstream>
#include <iostream>
#include <map>
#include <string>

#include "MRCPP/mrcpp_declarations.h"

namespace mrcpp {

/** @class Plotter
 *
 * @brief Class for plotting multivariate functions
 *
 * This class will generate an equidistant grid in one (line), two (surf)
 * or three (cube) dimensions, and subsequently evaluate the function on
 * this grid.
 *
 * The grid is generated from the vectors A, B and C, relative to the origin O:
 *  - a linePlot will plot the line spanned by A, starting from O
 *  - a surfPlot will plot the area spanned by A and B, starting from O
 *  - a cubePlot will plot the volume spanned by A, B and C, starting from O
 *
 * The vectors A, B and C do not necessarily have to be orthogonal.
 *
 * The parameter `D` refers to the dimension of the _function_, not the
 * dimension of the plot.
 *
 */

template <int D> class Plotter {
public:
    explicit Plotter(const Coord<D> &o = {});
    virtual ~Plotter() = default;

    void setSuffix(int t, const std::string &s);
    void setOrigin(const Coord<D> &o);
    void setRange(const Coord<D> &a, const Coord<D> &b = {}, const Coord<D> &c = {});

    void gridPlot(const MWTree<D> &tree, const std::string &fname);
    void linePlot(const std::array<int, 1> &npts, const RepresentableFunction<D> &func, const std::string &fname);
    void surfPlot(const std::array<int, 2> &npts, const RepresentableFunction<D> &func, const std::string &fname);
    void cubePlot(const std::array<int, 3> &npts, const RepresentableFunction<D> &func, const std::string &fname);

    enum type { Line, Surface, Cube, Grid };

protected:
    Coord<D> O{}; // Plot origin
    Coord<D> A{}; // Vector for line plot
    Coord<D> B{}; // Vector for surf plot
    Coord<D> C{}; // Vector for cube plot
    std::ofstream fstrm{};
    std::ofstream *fout{nullptr};
    std::map<int, std::string> suffix{};

    Coord<D> calcStep(const Coord<D> &vec, int pts) const;
    Eigen::MatrixXd calcLineCoordinates(int pts_a) const;
    Eigen::MatrixXd calcSurfCoordinates(int pts_a, int pts_b) const;
    Eigen::MatrixXd calcCubeCoordinates(int pts_a, int pts_b, int pts_c) const;

    Eigen::VectorXd evaluateFunction(const RepresentableFunction<D> &func, const Eigen::MatrixXd &coords) const;

    void writeData(const Eigen::MatrixXd &coords, const Eigen::VectorXd &values);
    virtual void writeCube(const std::array<int, 3> &npts, const Eigen::VectorXd &values);

    void writeGrid(const MWTree<D> &tree);
    void writeNodeGrid(const MWNode<D> &node, const std::string &color);

private:
    bool verifyRange(int dim) const;
    void openPlot(const std::string &fname);
    void closePlot();
};

} // namespace mrcpp
