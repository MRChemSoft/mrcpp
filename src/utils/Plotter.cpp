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

#include "Plotter.h"
#include "Printer.h"
#include "functions/RepresentableFunction.h"
#include "math_utils.h"
#include "trees/MWNode.h"

using namespace Eigen;

namespace mrcpp {

template <int D, typename T>
Plotter<D, T>::Plotter(const Coord<D> &o)
        : O(o) {
    setSuffix(Plotter<D, T>::Line, ".line");
    setSuffix(Plotter<D, T>::Surface, ".surf");
    setSuffix(Plotter<D, T>::Cube, ".cube");
    setSuffix(Plotter<D, T>::Grid, ".grid");
}

template <int D, typename T> void Plotter<D, T>::setSuffix(int t, const std::string &s) {
    this->suffix.insert(std::pair<int, std::string>(t, s));
}

template <int D, typename T> void Plotter<D, T>::setOrigin(const Coord<D> &o) {
    this->O = o;
}

template <int D, typename T> void Plotter<D, T>::setRange(const Coord<D> &a, const Coord<D> &b, const Coord<D> &c) {
    this->A = a;
    this->B = b;
    this->C = c;
}

template <int D, typename T> void Plotter<D, T>::gridPlot(const MWTree<D, T> &tree, const std::string &fname) {
    println(20, "----------Grid Plot-----------");
    std::stringstream file;
    file << fname << this->suffix[Plotter<D, T>::Grid];
    openPlot(file.str());
    writeGrid(tree);
    closePlot();
    printout(20, std::endl);
}

template <int D, typename T> void Plotter<D, T>::linePlot(const std::array<int, 1> &npts, const RepresentableFunction<D, T> &func, const std::string &fname) {
    println(20, "----------Line Plot-----------");
    std::stringstream file;
    file << fname << this->suffix[Plotter<D, T>::Line];
    if (verifyRange(1)) {
        Eigen::MatrixXd coords = calcLineCoordinates(npts[0]);
        Eigen::Matrix<T, Eigen::Dynamic, 1> values = evaluateFunction(func, coords);
        openPlot(file.str());
        writeData(coords, values);
        closePlot();
    } else {
        MSG_ERROR("Zero range");
    }
    printout(20, std::endl);
}

template <int D, typename T> void Plotter<D, T>::surfPlot(const std::array<int, 2> &npts, const RepresentableFunction<D, T> &func, const std::string &fname) {
    println(20, "--------Surface Plot----------");
    std::stringstream file;
    file << fname << this->suffix[Plotter<D, T>::Surface];
    if (verifyRange(2)) {
        Eigen::MatrixXd coords = calcSurfCoordinates(npts[0], npts[1]);
        Eigen::Matrix<T, Eigen::Dynamic, 1> values = evaluateFunction(func, coords);
        openPlot(file.str());
        writeData(coords, values);
        closePlot();
    } else {
        MSG_ERROR("Zero range");
    }
    printout(20, std::endl);
}

template <int D, typename T> void Plotter<D, T>::cubePlot(const std::array<int, 3> &npts, const RepresentableFunction<D, T> &func, const std::string &fname) {
    println(20, "----------Cube Plot-----------");
    std::stringstream file;
    file << fname << this->suffix[Plotter<D, T>::Cube];
    if (verifyRange(3)) {
        Eigen::MatrixXd coords = calcCubeCoordinates(npts[0], npts[1], npts[2]);
        Eigen::Matrix<T, Eigen::Dynamic, 1> values = evaluateFunction(func, coords);
        openPlot(file.str());
        writeCube(npts, values);
        closePlot();
    } else {
        MSG_ERROR("Zero range");
    }
    printout(20, std::endl);
}

template <int D, typename T> Eigen::MatrixXd Plotter<D, T>::calcLineCoordinates(int pts_a) const {
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

template <int D, typename T> Eigen::MatrixXd Plotter<D, T>::calcSurfCoordinates(int pts_a, int pts_b) const {
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

template <int D, typename T> Eigen::MatrixXd Plotter<D, T>::calcCubeCoordinates(int pts_a, int pts_b, int pts_c) const {
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

template <int D, typename T> Eigen::Matrix<T, Eigen::Dynamic, 1> Plotter<D, T>::evaluateFunction(const RepresentableFunction<D, T> &func, const Eigen::MatrixXd &coords) const {
    auto npts = coords.rows();
    if (npts == 0) MSG_ERROR("Empty coordinates");
    Eigen::Matrix<T, Eigen::Dynamic, 1> values = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(npts);
#pragma omp parallel for schedule(static) num_threads(mrcpp_get_num_threads())
    for (auto i = 0; i < npts; i++) {
        Coord<D> r{};
        for (auto d = 0; d < D; d++) r[d] = coords(i, d);
        values[i] = func.evalf(r);
    }
    return values;
}

template <int D, typename T> void Plotter<D, T>::writeData(const Eigen::MatrixXd &coords, const Eigen::Matrix<T, Eigen::Dynamic, 1> &values) {
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

template <int D, typename T> void Plotter<D, T>::writeCube(const std::array<int, 3> &npts, const Eigen::Matrix<T, Eigen::Dynamic, 1> &values) {
    NOT_IMPLEMENTED_ABORT
}

template <int D, typename T> void Plotter<D, T>::writeNodeGrid(const MWNode<D, T> &node, const std::string &color) {
    NOT_IMPLEMENTED_ABORT
}

template <int D, typename T> void Plotter<D, T>::writeGrid(const MWTree<D, T> &tree) {
    NOT_IMPLEMENTED_ABORT
}

template <int D, typename T> void Plotter<D, T>::openPlot(const std::string &fname) {
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

template <int D, typename T> void Plotter<D, T>::closePlot() {
    if (this->fout != nullptr) this->fout->close();
    this->fout = nullptr;
}

template <> void Plotter<3>::writeCube(const std::array<int, 3> &npts, const Eigen::VectorXd &values) {
    std::ofstream &o = *this->fout;

    Coord<3> a = calcStep(this->A, npts[0]);
    Coord<3> b = calcStep(this->B, npts[1]);
    Coord<3> c = calcStep(this->C, npts[2]);

    o << "Cube file format" << std::endl;
    o << "Generated by MRCPP" << std::endl;

    o.setf(std::ios::scientific);
    o.precision(6);

    o << std::setw(5) << 0;
    o << std::setw(15) << this->O[0];
    o << std::setw(15) << this->O[1];
    o << std::setw(15) << this->O[2] << std::endl;

    o << std::setw(5) << npts[0];
    o << std::setw(15) << a[0];
    o << std::setw(15) << a[1];
    o << std::setw(15) << a[2] << std::endl;

    o << std::setw(5) << npts[1];
    o << std::setw(15) << b[0];
    o << std::setw(15) << b[1];
    o << std::setw(15) << b[2] << std::endl;

    o << std::setw(5) << npts[2];
    o << std::setw(15) << c[0];
    o << std::setw(15) << c[1];
    o << std::setw(15) << c[2] << std::endl;

    o.precision(4);
    for (int n = 0; n < values.size(); n++) {
        o << std::setw(12) << values[n];
        if (n % 6 == 5) o << std::endl;
    }
}

template <> void Plotter<3>::writeNodeGrid(const MWNode<3> &node, const std::string &color) {
    double origin[3] = {0, 0, 0};
    double length = std::pow(2.0, -node.getScale());
    std::ostream &o = *this->fout;

    for (int d = 0; d < 3; d++) origin[d] = node.getNodeIndex()[d] * length;

    o << origin[0] << " " << origin[1] << " " << origin[2] << " " << color << origin[0] << " " << origin[1] << " " << origin[2] + length << " " << color << origin[0] << " " << origin[1] + length
      << " " << origin[2] + length << " " << color << origin[0] << " " << origin[1] + length << " " << origin[2] << color << std::endl;

    o << origin[0] << " " << origin[1] << " " << origin[2] << " " << color << origin[0] << " " << origin[1] << " " << origin[2] + length << " " << color << origin[0] + length << " " << origin[1]
      << " " << origin[2] + length << " " << color << origin[0] + length << " " << origin[1] << " " << origin[2] << color << std::endl;
    o << origin[0] << " " << origin[1] << " " << origin[2] << " " << color << origin[0] << " " << origin[1] + length << " " << origin[2] << " " << color << origin[0] + length << " "
      << origin[1] + length << " " << origin[2] << " " << color << origin[0] + length << " " << origin[1] << " " << origin[2] << color << std::endl;

    o << origin[0] + length << " " << origin[1] + length << " " << origin[2] + length << " " << color << origin[0] + length << " " << origin[1] + length << " " << origin[2] << " " << color
      << origin[0] + length << " " << origin[1] << " " << origin[2] << " " << color << origin[0] + length << " " << origin[1] << " " << origin[2] + length << color << std::endl;

    o << origin[0] + length << " " << origin[1] + length << " " << origin[2] + length << " " << color << origin[0] + length << " " << origin[1] + length << " " << origin[2] << " " << color
      << origin[0] << " " << origin[1] + length << " " << origin[2] << " " << color << origin[0] << " " << origin[1] + length << " " << origin[2] + length << color << std::endl;

    o << origin[0] + length << " " << origin[1] + length << " " << origin[2] + length << " " << color << origin[0] + length << " " << origin[1] << " " << origin[2] + length << " " << color
      << origin[0] << " " << origin[1] << " " << origin[2] + length << " " << color << origin[0] << " " << origin[1] + length << " " << origin[2] + length << color << std::endl;
}

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

template <int D, typename T> bool Plotter<D, T>::verifyRange(int dim) const {
    auto is_len_zero = [](Coord<D> vec) {
        double vec_sq = 0.0;
        for (auto d = 0; d < D; d++) vec_sq += vec[d] * vec[d];
        if (std::sqrt(vec_sq) < MachineZero) return true;
        return false;
    };

    if (is_len_zero(this->A)) return false;
    if (dim == 2 or dim == 3) {
        if (is_len_zero(this->B)) return false;
    }
    if (dim == 3) {
        if (is_len_zero(this->C)) return false;
    }

    return true;
}

template <int D, typename T> Coord<D> Plotter<D, T>::calcStep(const Coord<D> &vec, int pts) const {
    Coord<D> step;
    for (auto d = 0; d < D; d++) step[d] = vec[d] / (pts - 1.0);
    return step;
}

template class Plotter<1, double>;
template class Plotter<2, double>;
template class Plotter<3, double>;

template class Plotter<1, ComplexDouble>;
template class Plotter<2, ComplexDouble>;
template class Plotter<3, ComplexDouble>;

} // namespace mrcpp
