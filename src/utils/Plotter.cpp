/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2020 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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

#include "Plotter.h"
#include "Printer.h"
#include "functions/RepresentableFunction.h"
#include "math_utils.h"
#include "trees/MWNode.h"

using namespace Eigen;

namespace mrcpp {

/** Plotter constructor

    Arguments:
    o:                 origin, default (0, 0, ... , 0)
*/
template <int D>
Plotter<D>::Plotter(const Coord<D> &o)
        : O(o) {
    setSuffix(Plotter<D>::Line, ".line");
    setSuffix(Plotter<D>::Surface, ".surf");
    setSuffix(Plotter<D>::Cube, ".cube");
    setSuffix(Plotter<D>::Grid, ".grid");
}

/** Set file extension for output file

    The file name you decide for the output will get a predefined suffix
    that differentiates between different types of plot.

    Default values:
    line:    ".line"
    surface: ".surf"
    cube:    ".cube"
    grid:    ".grid"
*/
template <int D> void Plotter<D>::setSuffix(int t, const std::string &s) {
    this->suffix.insert(std::pair<int, std::string>(t, s));
}

/** Grid plot of a MWTree

    Writes a file named fname + file extension (".grid" as default) to be read
    by geomview to visualize the grid (of endNodes)	where the multiresolution
    function is defined. In MPI, each process will write a separate file, and
    will print only	nodes owned by itself (pluss the rootNodes).
*/
template <int D> void Plotter<D>::gridPlot(const MWTree<D> &tree, const std::string &fname) {
    println(20, "----------Grid Plot-----------");
    std::stringstream file;
    file << fname << this->suffix[Plotter<D>::Grid];
    openPlot(file.str());
    writeGrid(tree);
    closePlot();
    printout(20, std::endl);
}

/** Parametric plot of a function

    Plots the function func parametrically along the vector A starting
    from the origin O to a file named fname + file extension (".line"
    as default).
*/
template <int D>
void Plotter<D>::linePlot(const std::array<int, 1> &npts,
                          const RepresentableFunction<D> &func,
                          const std::string &fname) {
    println(20, "----------Line Plot-----------");
    std::stringstream file;
    file << fname << this->suffix[Plotter<D>::Line];
    if (verifyRange(1)) { // Verifies only A vector
        Eigen::MatrixXd coords = calcLineCoordinates(npts[0]);
        Eigen::VectorXd values = evaluateFunction(func, coords);
        openPlot(file.str());
        writeData(coords, values);
        closePlot();
    } else {
        MSG_ERROR("Zero range");
    }
    printout(20, std::endl);
}

/** Surface plot of a function

    Plots the function func in 2D on the area spanned by the two vectors
    A and B, starting from the origin O, to a file named fname + file
    extension (".surf" as default).
*/
template <int D>
void Plotter<D>::surfPlot(const std::array<int, 2> &npts,
                          const RepresentableFunction<D> &func,
                          const std::string &fname) {
    println(20, "--------Surface Plot----------");
    std::stringstream file;
    file << fname << this->suffix[Plotter<D>::Surface];
    if (verifyRange(2)) { // Verifies A and B vectors
        Eigen::MatrixXd coords = calcSurfCoordinates(npts[0], npts[1]);
        Eigen::VectorXd values = evaluateFunction(func, coords);
        openPlot(file.str());
        writeData(coords, values);
        closePlot();
    } else {
        MSG_ERROR("Zero range");
    }
    printout(20, std::endl);
}

/** Cubic plot of a function

    Plots the function func in 3D in the volume spanned by the three vectors
    A, B and C, starting from the origin O, to a file named fname + file
    extension (".cube" as default).
*/
template <int D>
void Plotter<D>::cubePlot(const std::array<int, 3> &npts,
                          const RepresentableFunction<D> &func,
                          const std::string &fname) {
    println(20, "----------Cube Plot-----------");
    std::stringstream file;
    file << fname << this->suffix[Plotter<D>::Cube];
    if (verifyRange(3)) { // Verifies A, B and C vectors
        Eigen::MatrixXd coords = calcCubeCoordinates(npts[0], npts[1], npts[2]);
        Eigen::VectorXd values = evaluateFunction(func, coords);
        openPlot(file.str());
        writeCube(npts, values);
        closePlot();
    } else {
        MSG_ERROR("Zero range");
    }
    printout(20, std::endl);
}

/** Calculating coordinates to be evaluated

    Generating a vector of pts_a equidistant coordinates that makes up the
    vector A in D dimensions, starting from the origin O.
*/
template <int D> Eigen::MatrixXd Plotter<D>::calcLineCoordinates(int pts_a) const {
    MatrixXd coords;
    if (pts_a > 0) {
        Coord<D> a = calcStep(this->A, pts_a);
        coords = MatrixXd::Zero(pts_a, D);
        for (auto i = 0; i < pts_a; i++) {
            for (auto d = 0; d < D; d++) coords(i, d) = this->O[d] + i * a[d];
        }
    } else {
        MSG_ERROR("Invalid number of points for plotting");
    }
    return coords;
}

/** Calculating coordinates to be evaluated

    Generating a vector of equidistant coordinates that makes up the
    area spanned by vectors A and B in D dimensions, starting from the
    origin O.
*/
template <int D> Eigen::MatrixXd Plotter<D>::calcSurfCoordinates(int pts_a, int pts_b) const {
    if (D < 2) MSG_ERROR("Cannot surfPlot less than 2D");

    MatrixXd coords;
    int npts = pts_a * pts_b;
    if (npts > 0) {
        Coord<D> a = calcStep(this->A, pts_a);
        Coord<D> b = calcStep(this->B, pts_b);

        auto n = 0;
        coords = MatrixXd::Zero(npts, D);
        for (auto i = 0; i < pts_a; i++) {
            for (auto j = 0; j < pts_b; j++) {
                for (auto d = 0; d < D; d++) coords(n, d) = this->O[d] + i * a[d] + j * b[d];
                n++;
            }
        }
    } else {
        MSG_ERROR("No points to plot");
    }
    return coords;
}

/** Calculating coordinates to be evaluated

    Generating a vector of equidistant coordinates that makes up the
    volume spanned by vectors A, B and C in D dimensions, starting from
    the origin O.
*/
template <int D> Eigen::MatrixXd Plotter<D>::calcCubeCoordinates(int pts_a, int pts_b, int pts_c) const {
    if (D < 3) MSG_ERROR("Cannot cubePlot less than 3D function");

    MatrixXd coords;
    int npts = pts_a * pts_b * pts_c;
    if (npts > 0) {
        Coord<D> a = calcStep(this->A, pts_a);
        Coord<D> b = calcStep(this->B, pts_b);
        Coord<D> c = calcStep(this->C, pts_c);

        auto n = 0;
        coords = MatrixXd::Zero(npts, D);
        for (auto i = 0; i < pts_a; i++) {
            for (auto j = 0; j < pts_b; j++) {
                for (auto k = 0; k < pts_c; k++) {
                    for (auto d = 0; d < D; d++) coords(n, d) = this->O[d] + i * a[d] + j * b[d] + k * c[d];
                    n++;
                }
            }
        }
    } else {
        MSG_ERROR("No points to plot");
    }
    return coords;
}

/** Evaluating a function in a set of predfined coordinates

    Given that the set of coordinates ("coords") has been calculated, this
    routine evaluates the function in these points and stores the results
    in the vector "values".
*/
template <int D>
Eigen::VectorXd Plotter<D>::evaluateFunction(const RepresentableFunction<D> &func,
                                             const Eigen::MatrixXd &coords) const {
    auto npts = coords.rows();
    if (npts == 0) MSG_ERROR("Empty coordinates");
    Eigen::VectorXd values = VectorXd::Zero(npts);
#pragma omp parallel for schedule(static)
    for (auto i = 0; i < npts; i++) {
        Coord<D> r{};
        for (auto d = 0; d < D; d++) r[d] = coords(i, d);
        values[i] = func.evalf(r);
    }
    return values;
}

/** Writing plot data to file

    This will write the contents of the "coords" matrix along with the function
    values to the file stream fout.	File will contain on each line the point
    number (between 0 and nPoints),	coordinates 1 through D and the function
    value.
*/
template <int D> void Plotter<D>::writeData(const Eigen::MatrixXd &coords, const Eigen::VectorXd &values) {
    if (coords.rows() != values.size()) INVALID_ARG_ABORT;
    std::ofstream &o = *this->fout;
    for (auto i = 0; i < values.size(); i++) {
        o.precision(8);
        o.setf(std::ios::showpoint);
        for (auto d = 0; d < D; d++) o << coords(i, d) << " ";
        o.precision(12);
        o << values[i];
        o << std::endl;
    }
}

// Specialized for D=3 below
template <int D> void Plotter<D>::writeCube(const std::array<int, 3> &npts, const Eigen::VectorXd &values) {
    NOT_IMPLEMENTED_ABORT
}

// Specialized for D=3 below
template <int D> void Plotter<D>::writeNodeGrid(const MWNode<D> &node, const std::string &color) {
    NOT_IMPLEMENTED_ABORT
}

// Specialized for D=3 below
template <int D> void Plotter<D>::writeGrid(const MWTree<D> &tree) {
    NOT_IMPLEMENTED_ABORT
}

/** Opening file for output

    Opens a file output stream fout for file named fname.
*/
template <int D> void Plotter<D>::openPlot(const std::string &fname) {
    if (fname.empty()) {
        if (this->fout == nullptr) {
            MSG_ERROR("Plot file not set!");
            return;
        } else if (this->fout->fail()) {
            MSG_ERROR("Plot file not set!");
            return;
        }
    } else {
        if (this->fout != nullptr) this->fout->close();
        this->fout = &this->fstrm;
        this->fout->open(fname.c_str());
        if (this->fout->bad()) {
            MSG_ERROR("File error");
            return;
        }
    }
}

/** Closing file

    Closes the file output stream fout.
*/
template <int D> void Plotter<D>::closePlot() {
    if (this->fout != nullptr) this->fout->close();
    this->fout = nullptr;
}

/** Writing plot data to file

    This will write a cube file (readable by blob) of the function values
    previously calculated (the "values" vector).
*/
template <> void Plotter<3>::writeCube(const std::array<int, 3> &npts, const Eigen::VectorXd &values) {
    std::ofstream &o = *this->fout;

    Coord<3> a = calcStep(this->A, npts[0]);
    Coord<3> b = calcStep(this->B, npts[1]);
    Coord<3> c = calcStep(this->C, npts[2]);

    o << "Cube file format" << std::endl;
    o << "Generated by MRCPP" << std::endl;

    o.setf(std::ios::scientific);
    o.precision(6);

    // Origin
    o << std::setw(5) << 0;
    o << std::setw(15) << this->O[0];
    o << std::setw(15) << this->O[1];
    o << std::setw(15) << this->O[2] << std::endl;

    // Vector A
    o << std::setw(5) << npts[0];
    o << std::setw(15) << a[0];
    o << std::setw(15) << a[1];
    o << std::setw(15) << a[2] << std::endl;

    // Vector B
    o << std::setw(5) << npts[1];
    o << std::setw(15) << b[0];
    o << std::setw(15) << b[1];
    o << std::setw(15) << b[2] << std::endl;

    // Vector C
    o << std::setw(5) << npts[2];
    o << std::setw(15) << c[0];
    o << std::setw(15) << c[1];
    o << std::setw(15) << c[2] << std::endl;

    // Function values
    o.precision(4);
    for (int n = 0; n < values.size(); n++) {
        o << std::setw(12) << values[n];
        if (n % 6 == 5) o << std::endl; // Line break after 6 values
    }
}

template <> void Plotter<3>::writeNodeGrid(const MWNode<3> &node, const std::string &color) {
    double origin[3] = {0, 0, 0};
    double length = std::pow(2.0, -node.getScale());
    std::ostream &o = *this->fout;

    for (int d = 0; d < 3; d++) origin[d] = node.getTranslation()[d] * length;

    o << origin[0] << " " << origin[1] << " " << origin[2] << " " << color << origin[0] << " " << origin[1] << " "
      << origin[2] + length << " " << color << origin[0] << " " << origin[1] + length << " " << origin[2] + length
      << " " << color << origin[0] << " " << origin[1] + length << " " << origin[2] << color << std::endl;

    o << origin[0] << " " << origin[1] << " " << origin[2] << " " << color << origin[0] << " " << origin[1] << " "
      << origin[2] + length << " " << color << origin[0] + length << " " << origin[1] << " " << origin[2] + length
      << " " << color << origin[0] + length << " " << origin[1] << " " << origin[2] << color << std::endl;
    o << origin[0] << " " << origin[1] << " " << origin[2] << " " << color << origin[0] << " " << origin[1] + length
      << " " << origin[2] << " " << color << origin[0] + length << " " << origin[1] + length << " " << origin[2] << " "
      << color << origin[0] + length << " " << origin[1] << " " << origin[2] << color << std::endl;

    o << origin[0] + length << " " << origin[1] + length << " " << origin[2] + length << " " << color
      << origin[0] + length << " " << origin[1] + length << " " << origin[2] << " " << color << origin[0] + length
      << " " << origin[1] << " " << origin[2] << " " << color << origin[0] + length << " " << origin[1] << " "
      << origin[2] + length << color << std::endl;

    o << origin[0] + length << " " << origin[1] + length << " " << origin[2] + length << " " << color
      << origin[0] + length << " " << origin[1] + length << " " << origin[2] << " " << color << origin[0] << " "
      << origin[1] + length << " " << origin[2] << " " << color << origin[0] << " " << origin[1] + length << " "
      << origin[2] + length << color << std::endl;

    o << origin[0] + length << " " << origin[1] + length << " " << origin[2] + length << " " << color
      << origin[0] + length << " " << origin[1] << " " << origin[2] + length << " " << color << origin[0] << " "
      << origin[1] << " " << origin[2] + length << " " << color << origin[0] << " " << origin[1] + length << " "
      << origin[2] + length << color << std::endl;
}

/** Writing grid data to file

    This will write a grid file (readable by geomview) of the grid (of endNodes)
    where the multiresolution function is defined.
    Currently only working in 3D.
*/
template <> void Plotter<3>::writeGrid(const MWTree<3> &tree) {
    std::ostream &o = *this->fout;
    o << "CQUAD" << std::endl;
    o.precision(6);
    std::string rootColor = " 1 1 1 0 ";
    std::string color = " 0 0 1 1 ";
    for (auto i = 0; i < tree.getRootBox().size(); i++) {
        const MWNode<3> &rootNode = tree.getRootMWNode(i);
        writeNodeGrid(rootNode, rootColor);
    }
    for (auto i = 0; i < tree.getNEndNodes(); i++) {
        const MWNode<3> &node = tree.getEndMWNode(i);
        writeNodeGrid(node, color);
    }
}

/** Checks the validity of the plotting range
 */
template <int D> bool Plotter<D>::verifyRange(int dim) const {
    if (dim == 1) {
        double len_a = std::sqrt(this->A[0] * this->A[0]);
        if (len_a < MachineZero) return false;
    }
    if (dim == 2) {
        double len_a = std::sqrt(A[0] * A[0] + A[1] * A[1]);
        double len_b = std::sqrt(B[0] * B[0] + B[1] * B[1]);
        if (len_a < MachineZero and len_b < MachineZero) return false;
    }
    if (dim == 3) {
        double len_a = std::sqrt(A[0] * A[0] + A[1] * A[1] + A[2] * A[2]);
        double len_b = std::sqrt(B[0] * B[0] + B[1] * B[1] + B[2] * B[2]);
        double len_c = std::sqrt(C[0] * C[0] + C[1] * C[1] + C[2] * C[2]);
        if ((len_a < MachineZero) and (len_b < MachineZero) and (len_c < MachineZero)) return false;
    }
    return true;
}

template <int D> Coord<D> Plotter<D>::calcStep(const Coord<D> &vec, int pts) const {
    Coord<D> step;
    for (auto d = 0; d < D; d++) step[d] = vec[d] / (pts - 1.0);
    return step;
}

template class Plotter<1>;
template class Plotter<2>;
template class Plotter<3>;

} // namespace mrcpp
