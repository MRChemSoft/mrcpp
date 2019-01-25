/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2019 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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

#pragma GCC system_header
#include <Eigen/Core>

#include <fstream>
#include <iostream>
#include <map>
#include <string>

#include "mrcpp_declarations.h"

namespace mrcpp {

template <int D> class Plotter {
public:
    Plotter(int npts = 1000, const double *a = 0, const double *b = 0);
    virtual ~Plotter() = default;

    void setRange(const double *a, const double *b);
    void setNPoints(int npts);
    void setSuffix(int t, const std::string &s);

    void linePlot(const RepresentableFunction<D> &func, const std::string &fname);
    void surfPlot(const RepresentableFunction<D> &func, const std::string &fname);
    void cubePlot(const RepresentableFunction<D> &func, const std::string &fname);

    void linePlot(FunctionTree<D> &func, const std::string &fname);
    void surfPlot(FunctionTree<D> &func, const std::string &fname);
    void cubePlot(FunctionTree<D> &func, const std::string &fname);
    void gridPlot(const MWTree<D> &tree, const std::string &fname);

    Eigen::VectorXd &linePlot(RepresentableFunction<D> &func);
    Eigen::VectorXd &surfPlot(RepresentableFunction<D> &func);
    Eigen::VectorXd &cubePlot(RepresentableFunction<D> &func);

    enum type { Line, Surface, Cube, Grid };

protected:
    std::ofstream fstrm;
    std::ofstream *fout;
    int nPoints;
    double A[D]; ///< lower left corner
    double B[D]; ///< upper right corner
    std::map<int, std::string> suffix;
    Eigen::MatrixXd coords;
    Eigen::VectorXd values;

    void calcLineCoordinates();
    void calcSurfCoordinates();
    void calcCubeCoordinates();

    void evaluateFunction(const RepresentableFunction<D> &func);
    void evaluateFunction(FunctionTree<D> &tree);
    bool verifyRange();

    void writeLineData();
    void writeSurfData();
    virtual void writeCubeData();

    void writeGrid(const MWTree<D> &tree);
    void writeNodeGrid(const MWNode<D> &node, const std::string &color);

private:
    void openPlot(const std::string &fname);
    void closePlot();
};

} // namespace mrcpp
