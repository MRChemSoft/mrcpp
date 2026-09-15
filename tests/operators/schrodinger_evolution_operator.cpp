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

#include "catch2/catch_all.hpp"

#include "factory_functions.h"

#include <cmath>
#include <memory>

#include "functions/special_functions.h"
#include "operators/MWOperator.h"
#include "operators/TimeEvolutionOperator.h"
#include "treebuilders/add.h"
#include "treebuilders/apply.h"
#include "treebuilders/project.h"

namespace schrodinger_evolution_operator {

TEST_CASE("Apply Schrodinger's evolution operator", "[apply_schrodinger_evolution], [schrodinger_evolution_operator], [mw_operator]") {
    const auto min_scale = 0;
    const auto max_depth = 25;
    const auto order = 4;
    const auto prec = 1.0e-7;

    double t1 = 0.001;
    double delta_t = 0.03;
    double t2 = delta_t + t1;

    auto basis = mrcpp::LegendreBasis(order);
    auto world = mrcpp::BoundingBox<1>(min_scale);
    auto MRA = mrcpp::MultiResolutionAnalysis<1>(world, basis, max_depth);

    // One operator holding the whole semigroup as a complex operator tree
    mrcpp::TimeEvolutionOperator<1> Exp(MRA, prec, delta_t);
    REQUIRE(Exp.iscomplex());

    double sigma = 0.001;
    double x0 = 0.5;

    auto f = [sigma, x0, t = t1](const mrcpp::Coord<1> &r) -> ComplexDouble { return mrcpp::free_particle_analytical_solution(r[0], x0, t, sigma); };
    auto g = [sigma, x0, t = t2](const mrcpp::Coord<1> &r) -> ComplexDouble { return mrcpp::free_particle_analytical_solution(r[0], x0, t, sigma); };

    mrcpp::FunctionTree<1, ComplexDouble> f_tree(MRA);
    mrcpp::project<1, ComplexDouble>(prec, f_tree, f);
    mrcpp::FunctionTree<1, ComplexDouble> g_tree(MRA);
    mrcpp::project<1, ComplexDouble>(prec, g_tree, g);

    mrcpp::FunctionTree<1, ComplexDouble> fout_tree(MRA);
    mrcpp::apply<1, ComplexDouble>(prec, fout_tree, Exp, f_tree, -1, false);

    // Agreement with the analytic solution at the later time
    mrcpp::FunctionTree<1, ComplexDouble> error(MRA);
    mrcpp::add<1, ComplexDouble>(prec, error, {1.0, 0.0}, fout_tree, {-1.0, 0.0}, g_tree, -1, false, false);
    double tolerance = prec * prec / 25.0;
    REQUIRE(error.getSquareNorm() == Catch::Approx(0.0).margin(tolerance));

    // This propagator generates an imaginary component from a real state.
    auto real_in = [sigma, x0](const mrcpp::Coord<1> &r) -> ComplexDouble {
        double dx = r[0] - x0;
        return {std::exp(-dx * dx / (2.0 * sigma)), 0.0};
    };
    mrcpp::FunctionTree<1, ComplexDouble> real_tree(MRA);
    mrcpp::project<1, ComplexDouble>(prec, real_tree, real_in);

    // The caller owns the trees returned by Real() and Imag()
    std::unique_ptr<mrcpp::FunctionTree<1, double>> in_imag(real_tree.Imag());
    REQUIRE(in_imag->getSquareNorm() == Catch::Approx(0.0).margin(1.0e-20));

    mrcpp::FunctionTree<1, ComplexDouble> real_out(MRA);
    mrcpp::apply<1, ComplexDouble>(prec, real_out, Exp, real_tree, -1, false);

    std::unique_ptr<mrcpp::FunctionTree<1, double>> out_imag(real_out.Imag());
    // Above the projection-noise level for this prec
    REQUIRE(out_imag->getSquareNorm() > 1.0e-12);
}

TEST_CASE("Apply Schrodinger's evolution operator on a uniform grid", "[apply_schrodinger_uniform], [schrodinger_evolution_operator], [mw_operator]") {
    const auto min_scale = 0;
    const auto max_depth = 25;
    const auto order = 4;
    const auto prec = 1.0e-7;

    double t1 = 0.001;
    double delta_t = 0.03;
    double t2 = delta_t + t1;

    auto basis = mrcpp::LegendreBasis(order);
    auto world = mrcpp::BoundingBox<1>(min_scale);
    auto MRA = mrcpp::MultiResolutionAnalysis<1>(world, basis, max_depth);

    // The fixed-scale overload, bypassing the adaptive build
    mrcpp::TimeEvolutionOperator<1> Exp(MRA, prec, delta_t, 7);
    REQUIRE(Exp.iscomplex());

    double sigma = 0.001;
    double x0 = 0.5;

    auto f = [sigma, x0, t = t1](const mrcpp::Coord<1> &r) -> ComplexDouble { return mrcpp::free_particle_analytical_solution(r[0], x0, t, sigma); };
    auto g = [sigma, x0, t = t2](const mrcpp::Coord<1> &r) -> ComplexDouble { return mrcpp::free_particle_analytical_solution(r[0], x0, t, sigma); };

    mrcpp::FunctionTree<1, ComplexDouble> f_tree(MRA);
    mrcpp::project<1, ComplexDouble>(prec, f_tree, f);
    mrcpp::FunctionTree<1, ComplexDouble> g_tree(MRA);
    mrcpp::project<1, ComplexDouble>(prec, g_tree, g);

    mrcpp::FunctionTree<1, ComplexDouble> fout_tree(MRA);
    mrcpp::apply<1, ComplexDouble>(prec, fout_tree, Exp, f_tree, -1, false);

    mrcpp::FunctionTree<1, ComplexDouble> error(MRA);
    mrcpp::add<1, ComplexDouble>(prec, error, {1.0, 0.0}, fout_tree, {-1.0, 0.0}, g_tree, -1, false, false);
    double tolerance = prec * prec / 25.0;
    REQUIRE(error.getSquareNorm() == Catch::Approx(0.0).margin(tolerance));
}

} // namespace schrodinger_evolution_operator
