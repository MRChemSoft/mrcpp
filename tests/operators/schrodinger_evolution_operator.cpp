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

#include "functions/GaussFunc.h"
#include "functions/special_functions.h"
#include "operators/MWOperator.h"
#include "operators/TimeEvolutionOperator.h"
#include "treebuilders/add.h"
#include "treebuilders/complex_apply.h"
#include "treebuilders/project.h"

namespace schrodinger_evolution_operator {

TEST_CASE("Apply Schrodinger's evolution operator", "[apply_schrodinger_evolution], [schrodinger_evolution_operator], [mw_operator]") {
    const auto min_scale = 0;
    const auto max_depth = 25;

    const auto order = 4;
    const auto prec = 1.0e-7;

    int finest_scale = 7; // for time evolution operator construction (not recommended to use more than 10)
    // int max_Jpower = 20;  //the amount of J integrals to be used in construction (20 should be enough)

    // Time moments:
    double t1 = 0.001;        // initial time moment (not recommended to use more than 0.001)
    double delta_t = 0.03;    // time step (not recommended to use less than 0.001)
    double t2 = delta_t + t1; // final time moment

    // Initialize world in the unit cube [0,1]
    auto basis = mrcpp::LegendreBasis(order);
    auto world = mrcpp::BoundingBox<1>(min_scale);
    auto MRA = mrcpp::MultiResolutionAnalysis<1>(world, basis, max_depth);

    // Time evolution operatror Exp(delta_t)
    mrcpp::TimeEvolutionOperator<1> ReExp(MRA, prec, delta_t, finest_scale, false);
    mrcpp::TimeEvolutionOperator<1> ImExp(MRA, prec, delta_t, finest_scale, true);

    // Analytical solution parameters for psi(x, t)
    double sigma = 0.001;
    double x0 = 0.5;

    // Functions f(x) = psi(x, t1) and g(x) = psi(x, t2)
    auto Re_f = [sigma, x0, t = t1](const mrcpp::Coord<1> &r) -> double { return mrcpp::free_particle_analytical_solution(r[0], x0, t, sigma).real(); };
    auto Im_f = [sigma, x0, t = t1](const mrcpp::Coord<1> &r) -> double { return mrcpp::free_particle_analytical_solution(r[0], x0, t, sigma).imag(); };
    auto Re_g = [sigma, x0, t = t2](const mrcpp::Coord<1> &r) -> double { return mrcpp::free_particle_analytical_solution(r[0], x0, t, sigma).real(); };
    auto Im_g = [sigma, x0, t = t2](const mrcpp::Coord<1> &r) -> double { return mrcpp::free_particle_analytical_solution(r[0], x0, t, sigma).imag(); };

    // Projecting functions
    mrcpp::FunctionTree<1> Re_f_tree(MRA);
    mrcpp::project<1, double>(prec, Re_f_tree, Re_f);
    mrcpp::FunctionTree<1> Im_f_tree(MRA);
    mrcpp::project<1, double>(prec, Im_f_tree, Im_f);
    mrcpp::FunctionTree<1> Re_g_tree(MRA);
    mrcpp::project<1, double>(prec, Re_g_tree, Re_g);
    mrcpp::FunctionTree<1> Im_g_tree(MRA);
    mrcpp::project<1, double>(prec, Im_g_tree, Im_g);

    // Output function trees
    mrcpp::FunctionTree<1> Re_fout_tree(MRA);
    mrcpp::FunctionTree<1> Im_fout_tree(MRA);

    // Complex objects for use in apply()
    mrcpp::ComplexObject<mrcpp::ConvolutionOperator<1>> E(ReExp, ImExp);
    mrcpp::ComplexObject<mrcpp::FunctionTree<1>> input(Re_f_tree, Im_f_tree);
    mrcpp::ComplexObject<mrcpp::FunctionTree<1>> output(Re_fout_tree, Im_fout_tree);

    // Apply operator Exp(delta_t) f(x)
    mrcpp::apply(prec, output, E, input);

    // Check g(x) = Exp(delta_t) f(x)
    mrcpp::FunctionTree<1> Re_error(MRA); // = Re_fout_tree - Re_g_tree
    mrcpp::FunctionTree<1> Im_error(MRA); // = Im_fout_tree - Im_g_tree

    // Re_error = Re_fout_tree - Re_g_tree
    mrcpp::add(prec, Re_error, 1.0, Re_fout_tree, -1.0, Re_g_tree);
    auto Re_sq_norm = Re_error.getSquareNorm(); // 1.7e-16

    // Im_error = Im_fout_tree - Im_g_tree
    mrcpp::add(prec, Im_error, 1.0, Im_fout_tree, -1.0, Im_g_tree);
    auto Im_sq_norm = Im_error.getSquareNorm(); // 1.7e-17

    double tolerance = prec * prec / 50.0; // 2.0e-16

    // std::cout << "Re_sq_norm = " << Re_sq_norm << std::endl;
    // std::cout << "Im_sq_norm = " << Im_sq_norm << std::endl;
    // std::cout << "tolerance = " << tolerance << std::endl;

    REQUIRE(Re_sq_norm == Catch::Approx(0.0).margin(tolerance));
    REQUIRE(Im_sq_norm == Catch::Approx(0.0).margin(tolerance));
}

} // namespace schrodinger_evolution_operator
