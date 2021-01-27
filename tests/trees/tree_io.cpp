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

#include "catch.hpp"

#include "factory_functions.h"

#include "treebuilders/grid.h"
#include "treebuilders/project.h"

using namespace mrcpp;

namespace tree_io {

template <int D> void testProjectFunction();

SCENARIO("FunctionTree IO", "[tree_io], [trees]") {
    const double prec = 1.0e-4;

    GaussFunc<3> *func = nullptr;
    initialize(&func);
    MultiResolutionAnalysis<3> *mra = nullptr;
    initialize(&mra);

    FunctionTree<3> f_tree(*mra);
    project(prec, f_tree, *func);

    const double ref_charge = f_tree.integrate();
    const double ref_norm = f_tree.getSquareNorm();
    const double ref_nodes = f_tree.getNNodes();

    WHEN("a function is saved") {
        f_tree.saveTree("f");
        THEN("the old tree remains unchanged") {
            REQUIRE(f_tree.integrate() == Approx(ref_charge).epsilon(1.0e-12));
            REQUIRE(f_tree.getSquareNorm() == Approx(ref_norm).epsilon(1.0e-12));
            REQUIRE(f_tree.getNNodes() == ref_nodes);
        }
        AND_WHEN("the saved function is load into a new tree") {
            FunctionTree<3> g_tree(*mra);
            g_tree.loadTree("f");
            THEN("the new tree is identical to the old") {
                REQUIRE(g_tree.integrate() == Approx(ref_charge).epsilon(1.0e-12));
                REQUIRE(g_tree.getSquareNorm() == Approx(ref_norm).epsilon(1.0e-12));
                REQUIRE(g_tree.getNNodes() == ref_nodes);
            }
        }
    }
    WHEN("a function with unused memory chunks is saved") {
        f_tree.clear();                     // leaves memory chunks
        project(100 * prec, f_tree, *func); // should use less chunks

        const int refChunks = f_tree.getNChunks();
        const int refChunksUsed = f_tree.getNChunksUsed();
        REQUIRE(refChunksUsed < refChunks);

        f_tree.saveTree("f");
        THEN("the old tree remains unchanged") {
            int nChunks = f_tree.getProjectedNodeAllocator().getNChunks();
            int nChunksUsed = f_tree.getProjectedNodeAllocator().getNChunksUsed();
            REQUIRE(nChunks == refChunks);
            REQUIRE(nChunksUsed == refChunksUsed);
        }
        AND_WHEN("the saved function is load into a new tree") {
            FunctionTree<3> g_tree(*mra);
            g_tree.loadTree("f");
            THEN("the new tree has no empty chunks") {
                int nChunks = g_tree.getProjectedNodeAllocator().getNChunks();
                int nChunksUsed = g_tree.getProjectedNodeAllocator().getNChunksUsed();
                REQUIRE(nChunksUsed == nChunks);
                REQUIRE(nChunksUsed == refChunksUsed);
            }
        }
    }
    // Delete saved file
    remove("f.tree");

    finalize(&mra);
    finalize(&func);
}

} // namespace tree_io
