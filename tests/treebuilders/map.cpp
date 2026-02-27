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

#include "functions/GaussPoly.h"
#include "treebuilders/WaveletAdaptor.h"
#include "treebuilders/grid.h"
#include "treebuilders/map.h"
#include "treebuilders/project.h"

using namespace mrcpp;

namespace mapping {

template <int D> void testMapping();

SCENARIO("Map a MW tree", "[map], [tree_builder]") {
    GIVEN("One MW functions in 1D") {
        testMapping<1>();
    }
    GIVEN("One MW functions in 2D") {
        testMapping<2>();
    }
    GIVEN("One MW functions in 3D") {
        testMapping<3>();
    }
}

template <int D> void testMapping() {
    const double prec = 1.0e-4;
    double alpha = 1.0;
    double beta = 50.0;
    double beta_ref = 100.0;

    double pos_c[3] = {-0.25, 0.35, 1.05};
    auto pos = details::convert_to_std_array<double, D>(pos_c);

    GaussFunc<D> inp_func(beta, alpha, pos);
    GaussFunc<D> ref_func(beta_ref, alpha, pos);

    MultiResolutionAnalysis<D> *mra = nullptr;
    initialize<D>(&mra);

    // Initialize trees
    FunctionTree<D> inp_tree(*mra);
    FunctionTree<D> ref_tree(*mra);

    // Build empty grids
    build_grid(inp_tree, inp_func);
    build_grid(ref_tree, ref_func);

    // Project functions
    project(prec, inp_tree, inp_func);
    project(prec, ref_tree, ref_func);

    const double ref_int = ref_tree.integrate();
    const double ref_norm = ref_tree.getSquareNorm();
    const double inp_norm = inp_tree.getSquareNorm();

    FMap<double, double> fmap = [](double val) { return val * val; };

    WHEN("the function is mapped") {
        FunctionTree<D> out_tree(*mra);

        map(prec, out_tree, inp_tree, fmap);

        THEN("the MW map equals the analytic map") {
            double out_int = out_tree.integrate();
            double out_norm = out_tree.getSquareNorm();
            REQUIRE(out_int == Catch::Approx(ref_int));
            REQUIRE(out_norm == Catch::Approx(ref_norm));
            REQUIRE(inp_norm == Catch::Approx(out_int));
        }
    }
    WHEN("the functions is mapped in-place") {
        inp_tree.map(fmap);
        THEN("the first function equals the analytic product") {
            double inp_int = inp_tree.integrate();
            double inp_norm = inp_tree.getSquareNorm();
            REQUIRE(inp_int == Catch::Approx(ref_int));
            REQUIRE(inp_norm == Catch::Approx(ref_norm));
        }
    }
    finalize(&mra);
}

} // namespace mapping
