#include "catch.hpp"

#include "factory_functions.h"

#include "project.h"
#include "grid.h"
#include "SerialFunctionTree.h"

using namespace mrcpp;

namespace tree_io {

template<int D> void testProjectFunction();

SCENARIO("FunctionTree IO", "[tree_io], [trees]") {
    const double prec = 1.0e-4;

    GaussFunc<3> *func = 0;
    initialize(&func);
    MultiResolutionAnalysis<3> *mra = 0;
    initialize(&mra);

    FunctionTree<3> f_tree(*mra);
    project(prec, f_tree, *func);

    const double ref_charge = f_tree.integrate();
    const double ref_norm = f_tree.getSquareNorm();
    const double ref_nodes = f_tree.getNNodes();

    WHEN("a function is saved") {
        f_tree.saveTree("f");
        THEN("the old tree remains unchanged") {
            REQUIRE( (f_tree.integrate() == Approx(ref_charge).epsilon(1.0e-12)) );
            REQUIRE( (f_tree.getSquareNorm() == Approx(ref_norm).epsilon(1.0e-12)) );
            REQUIRE( (f_tree.getNNodes() == ref_nodes) );
        }
        AND_WHEN("the saved function is load into a new tree") {
            FunctionTree<3> g_tree(*mra);
            g_tree.loadTree("f");
            THEN("the new tree is identical to the old") {
                REQUIRE( (g_tree.integrate() == Approx(ref_charge).epsilon(1.0e-12)) );
                REQUIRE( (g_tree.getSquareNorm() == Approx(ref_norm).epsilon(1.0e-12)) );
                REQUIRE( (g_tree.getNNodes() == ref_nodes) );
            }
        }
    }
    WHEN("a function with unused memory chunks is saved") {
        f_tree.clear();                   // leaves memory chunks
        project(100*prec, f_tree, *func); // should use less chunks

        const int refChunks = f_tree.getSerialFunctionTree()->getNChunks();
        const int refChunksUsed = f_tree.getSerialFunctionTree()->getNChunksUsed();
        REQUIRE( (refChunksUsed < refChunks) );

        f_tree.saveTree("f");
        THEN("the old tree remains unchanged") {
            int nChunks = f_tree.getSerialFunctionTree()->getNChunks();
            int nChunksUsed = f_tree.getSerialFunctionTree()->getNChunksUsed();
            REQUIRE( (nChunks == refChunks) );
            REQUIRE( (nChunksUsed == refChunksUsed) );
        }
        AND_WHEN("the saved function is load into a new tree") {
            FunctionTree<3> g_tree(*mra);
            g_tree.loadTree("f");
            THEN("the new tree has no empty chunks") {
                int nChunks = g_tree.getSerialFunctionTree()->getNChunks();
                int nChunksUsed = g_tree.getSerialFunctionTree()->getNChunksUsed();
                REQUIRE( (nChunksUsed == nChunks) );
                REQUIRE( (nChunksUsed == refChunksUsed) );
            }
        }
    }
    // Delete saved file
    remove("f.tree");

    finalize(&mra);
    finalize(&func);
}

} // namespace
