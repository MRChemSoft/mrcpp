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
#include "operators/MWOperator.h"
#include "operators/PoissonKernel.h"
#include "operators/PoissonOperator.h"
#include "treebuilders/CrossCorrelationCalculator.h"
#include "treebuilders/OperatorAdaptor.h"
#include "treebuilders/TreeBuilder.h"
#include "treebuilders/apply.h"
#include "treebuilders/grid.h"
#include "treebuilders/multiply.h"
#include "treebuilders/project.h"
#include "trees/BandWidth.h"
#include "operators/HeatOperator.h"
#include "functions/special_functions.h"
#include "treebuilders/complex_apply.h"
#include "treebuilders/add.h"

//using namespace mrcpp;

namespace heat_evolution_operator {


TEST_CASE("Apply heat evolution operator", "[apply_heat_evolution], [heat_evolution_operator], [mw_operator]") {
    const auto min_scale = 0;
    const auto max_depth = 25;

    const auto order = 5;
    const auto prec = 1.0e-8;

    
    // Time moment:
    double delta_t = 0.0005;

    // Initialize world in the unit cube [0,1]
    auto basis = mrcpp::LegendreBasis(order);
    auto world = mrcpp::BoundingBox<1>(min_scale);
    auto MRA = mrcpp::MultiResolutionAnalysis<1>(world, basis, max_depth);

    // Time evolution operatror Exp(delta_t)
    mrcpp::HeatOperator<1> H(MRA, delta_t, prec);
    
    // Analytical solution parameters for psi(x, t)
    double sigma = 0.001;
    double x0 = 0.5;

    mrcpp::GaussFunc<1> f( 1.0 / sigma,  1.0, {x0});
    double beta = 1.0 / (4.0 * delta_t + sigma);
    mrcpp::GaussFunc<1> g( beta,  std::sqrt( sigma * beta ), {x0});

    // Projecting functions
    mrcpp::FunctionTree<1> f_tree(MRA);
    mrcpp::project<1>(prec, f_tree, f);
    mrcpp::FunctionTree<1> g_tree(MRA);
    mrcpp::project<1>(prec, g_tree, g);
    
    // Apply operator H = Exp(delta_t) f(x)
    mrcpp::FunctionTree<1> output(MRA);
    mrcpp::apply(prec, output, H, f_tree);

    // Check g(x) = Exp(delta_t) f(x)
    mrcpp::FunctionTree<1> error(MRA);
    mrcpp::add(prec, error, 1.0, output, -1.0, g_tree);
    auto sq_norm = error.getSquareNorm();    //1.6e-20

    double tolerance = prec * prec / 1000.0; //1.0e-19

    REQUIRE(sq_norm == Catch::Approx(0.0).margin(tolerance));
}


} // namespace schrodinger_evolution_operator