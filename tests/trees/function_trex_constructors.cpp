#include "catch2/catch_all.hpp"
#include "factory_functions.h"

using namespace mrcpp;

namespace function_tree_constructors {

    template <int D, typename T>
    void testBasicConstructor();

    template <int D, typename T>
    void testNamedConstructor();

    template <int D, typename T>
    void testSharedMemoryConstructor();

    template <int D>
    void testComplexConstructorFromTwoRealTrees();


    SCENARIO("FunctionTree basic constructor", "[function_tree_ctor_basic][function_trex_constructors][trees]") {
        GIVEN("a 1D real FunctionTree") { testBasicConstructor<1, double>(); }
        GIVEN("a 2D real FunctionTree") { testBasicConstructor<2, double>(); }
        GIVEN("a 3D real FunctionTree") { testBasicConstructor<3, double>(); }

        GIVEN("a 1D complex FunctionTree") { testBasicConstructor<1, ComplexDouble>(); }
        GIVEN("a 2D complex FunctionTree") { testBasicConstructor<2, ComplexDouble>(); }
        GIVEN("a 3D complex FunctionTree") { testBasicConstructor<3, ComplexDouble>(); }
        }

    SCENARIO("FunctionTree named constructor", "[function_tree_ctor_name][function_trex_constructors][trees]") {
        GIVEN("a 1D real FunctionTree with explicit name") { testNamedConstructor<1, double>(); }
        GIVEN("a 2D real FunctionTree with explicit name") { testNamedConstructor<2, double>(); }
        GIVEN("a 3D real FunctionTree with explicit name") { testNamedConstructor<3, double>(); }

        GIVEN("a 1D complex FunctionTree with explicit name") { testNamedConstructor<1, ComplexDouble>(); }
        GIVEN("a 2D complex FunctionTree with explicit name") { testNamedConstructor<2, ComplexDouble>(); }
        GIVEN("a 3D complex FunctionTree with explicit name") { testNamedConstructor<3, ComplexDouble>(); }
    }

    SCENARIO("FunctionTree shared-memory constructor with nullptr", "[function_trex_ctor_shmem][function_tree_constructors][trees]") {
        GIVEN("a 1D real FunctionTree with explicit nullptr shared memory") { testSharedMemoryConstructor<1, double>(); }
        GIVEN("a 2D real FunctionTree with explicit nullptr shared memory") { testSharedMemoryConstructor<2, double>(); }
        GIVEN("a 3D real FunctionTree with explicit nullptr shared memory") { testSharedMemoryConstructor<3, double>(); }

        GIVEN("a 1D complex FunctionTree with explicit nullptr shared memory") { testSharedMemoryConstructor<1, ComplexDouble>(); }
        GIVEN("a 2D complex FunctionTree with explicit nullptr shared memory") { testSharedMemoryConstructor<2, ComplexDouble>(); }
        GIVEN("a 3D complex FunctionTree with explicit nullptr shared memory") { testSharedMemoryConstructor<3, ComplexDouble>(); }
    }

    SCENARIO("Construct complex FunctionTree from two real FunctionTrees",
             "[function_tree_complex_constructor], [function_tree], [function_trex_constructors], [trees]") {
        GIVEN("two real trees in 1D") { testComplexConstructorFromTwoRealTrees<1>(); }
        GIVEN("two real trees in 2D") { testComplexConstructorFromTwoRealTrees<2>(); }
        GIVEN("two real trees in 3D") { testComplexConstructorFromTwoRealTrees<3>(); }
    }

    template <int D, typename T>
    void testBasicConstructor() {
        MultiResolutionAnalysis<D>* mra = nullptr;
        initialize(&mra);

        FunctionTree<D, T> tree(*mra);

        THEN("the tree starts without generated nodes") {
            REQUIRE(tree.getNGenNodes() == 0);
        }

        THEN("the generated-node allocator exists and starts empty") {
            REQUIRE(tree.getGenNodeAllocator().getNNodes() == 0);
        }

        THEN("the root box contains at least one root node") {
            REQUIRE(tree.getRootBox().size() > 0);
        }

        THEN("the tree starts undefined, with negative square norm") {
            REQUIRE(tree.getSquareNorm() < 0.0);
        }

        finalize(&mra);
    }

    template <int D, typename T>
    void testNamedConstructor() {
        MultiResolutionAnalysis<D>* mra = nullptr;
        initialize(&mra);

        const std::string name = "ctor_test";
        FunctionTree<D, T> tree(*mra, name);

        THEN("the tree starts without generated nodes") {
            REQUIRE(tree.getNGenNodes() == 0);
        }

        THEN("the tree name is stored") {
            REQUIRE(tree.getName() == name);
        }

        THEN("the tree starts undefined, with negative square norm") {
            REQUIRE(tree.getSquareNorm() < 0.0);
        }

        finalize(&mra);
    }

    template <int D, typename T>
    void testSharedMemoryConstructor() {
        MultiResolutionAnalysis<D>* mra = nullptr;
        initialize(&mra);

        SharedMemory<T>* sh_mem = nullptr;
        const std::string name = "ctor_shmem_test";
        FunctionTree<D, T> tree(*mra, sh_mem, name);

        THEN("the tree starts without generated nodes") {
            REQUIRE(tree.getNGenNodes() == 0);
        }

        THEN("the generated-node allocator exists and starts empty") {
            REQUIRE(tree.getGenNodeAllocator().getNNodes() == 0);
        }

        THEN("the tree name is stored") {
            REQUIRE(tree.getName() == name);
        }

        THEN("the tree starts undefined, with negative square norm") {
            REQUIRE(tree.getSquareNorm() < 0.0);
        }

        finalize(&mra);
    }

    template <int D>
    void testComplexConstructorFromTwoRealTrees() {
        MultiResolutionAnalysis<D>* mra = nullptr;
        initialize(&mra);

        FunctionTree<D, double> realTree(*mra);
        FunctionTree<D, double> imagTree(*mra);

        realTree.setZero();
        imagTree.setZero();

        FunctionTree<D, ComplexDouble> ctree(realTree, imagTree, nullptr, "complex_from_two_real");

        Coord<D> r;
        if (r.size() >= 1) r[0] = -0.20;
        if (r.size() >= 2) r[1] =  0.60;
       if (r.size() >= 3) r[2] =  0.76;

       THEN("the constructed complex tree evaluates to zero") {
           auto val = ctree.evalf(r);
            REQUIRE(val.real() == Catch::Approx(0.0));
            REQUIRE(val.imag() == Catch::Approx(0.0));
        }

        THEN("the constructed complex tree has zero squared norm") {
            REQUIRE(ctree.getSquareNorm() == Catch::Approx(0.0));
        }

        THEN("the constructed complex tree integrates to zero") {
            auto val = ctree.integrate();
            REQUIRE(val.real() == Catch::Approx(0.0));
            REQUIRE(val.imag() == Catch::Approx(0.0));
        }

        finalize(&mra);
    }
} // namespace function_tree_constructors