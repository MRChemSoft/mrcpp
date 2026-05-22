#include "catch2/catch_all.hpp"
#include "factory_functions.h"
#include "treebuilders/grid.h"
#include "treebuilders/project.h"

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

    template <int D>
    void makeNonZeroTree(FunctionTree<D, double>& tree,
                         double alpha,
                         double beta,
                         const Coord<D>& pos,
                         const std::array<int, static_cast<std::size_t>(D)>& power);

    template <int D>
    void checkComplexTreeMatchesParts(const FunctionTree<D, double>& realTree,
                                      const FunctionTree<D, double>& imagTree,
                                      const FunctionTree<D, ComplexDouble>& ctree);

    template <int D>
    void testComplexConstructor(bool useRealPart, bool useImagPart);

    SCENARIO("FunctionTree basic constructor", "[function_tree_constructors]") {
        GIVEN("a 1D real FunctionTree") { testBasicConstructor<1, double>(); }
        GIVEN("a 2D real FunctionTree") { testBasicConstructor<2, double>(); }
        GIVEN("a 3D real FunctionTree") { testBasicConstructor<3, double>(); }

        GIVEN("a 1D complex FunctionTree") { testBasicConstructor<1, ComplexDouble>(); }
        GIVEN("a 2D complex FunctionTree") { testBasicConstructor<2, ComplexDouble>(); }
        GIVEN("a 3D complex FunctionTree") { testBasicConstructor<3, ComplexDouble>(); }
        }

    SCENARIO("FunctionTree named constructor", "[function_tree_constructors]") {
        GIVEN("a 1D real FunctionTree with explicit name") { testNamedConstructor<1, double>(); }
        GIVEN("a 2D real FunctionTree with explicit name") { testNamedConstructor<2, double>(); }
        GIVEN("a 3D real FunctionTree with explicit name") { testNamedConstructor<3, double>(); }

        GIVEN("a 1D complex FunctionTree with explicit name") { testNamedConstructor<1, ComplexDouble>(); }
        GIVEN("a 2D complex FunctionTree with explicit name") { testNamedConstructor<2, ComplexDouble>(); }
        GIVEN("a 3D complex FunctionTree with explicit name") { testNamedConstructor<3, ComplexDouble>(); }
    }

    SCENARIO("FunctionTree shared-memory constructor with nullptr", "[function_tree_constructors]") {
        GIVEN("a 1D real FunctionTree with explicit nullptr shared memory") { testSharedMemoryConstructor<1, double>(); }
        GIVEN("a 2D real FunctionTree with explicit nullptr shared memory") { testSharedMemoryConstructor<2, double>(); }
        GIVEN("a 3D real FunctionTree with explicit nullptr shared memory") { testSharedMemoryConstructor<3, double>(); }

        GIVEN("a 1D complex FunctionTree with explicit nullptr shared memory") { testSharedMemoryConstructor<1, ComplexDouble>(); }
        GIVEN("a 2D complex FunctionTree with explicit nullptr shared memory") { testSharedMemoryConstructor<2, ComplexDouble>(); }
        GIVEN("a 3D complex FunctionTree with explicit nullptr shared memory") { testSharedMemoryConstructor<3, ComplexDouble>(); }
    }

    SCENARIO("Construct complex FunctionTree from two zero real FunctionTrees",
             "[function_tree_constructors]") {
        GIVEN("two real trees in 1D") { testComplexConstructorFromTwoRealTrees<1>(); }
        GIVEN("two real trees in 2D") { testComplexConstructorFromTwoRealTrees<2>(); }
        GIVEN("two real trees in 3D") { testComplexConstructorFromTwoRealTrees<3>(); }
    }

    SCENARIO("Construct complex FunctionTree from two real FunctionTrees with nonzero values",
             "[function_tree_constructors]") {
 
        GIVEN("real-only input") {
            testComplexConstructor<1>(true,  false);
            testComplexConstructor<2>(true,  false);
            testComplexConstructor<3>(true,  false);
        }

        GIVEN("imag-only input") {
            testComplexConstructor<1>(false, true);
            testComplexConstructor<2>(false, true);
            testComplexConstructor<3>(false, true);
        }

        GIVEN("both inputs nonzero") {
            testComplexConstructor<1>(true,  true);
            testComplexConstructor<2>(true,  true);
            testComplexConstructor<3>(true,  true);
        }
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
/*
        THEN("the tree starts undefined, with negative square norm") {
            REQUIRE(tree.getSquareNorm() < 0.0);
        }
*/
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
/*
        THEN("the tree starts undefined, with negative square norm") {
            REQUIRE(tree.getSquareNorm() < 0.0);
        }
*/
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

    template <int D>
    void makeNonZeroTree(FunctionTree<D, double>& tree,
                         double alpha,
                         double beta,
                         const Coord<D>& pos,
                         const std::array<int, static_cast<std::size_t>(D)>& power) {
        auto prec = 1.0e-4;
        GaussFunc<D> f(alpha, beta, pos, power);
        build_grid<D>(tree, f);
        project(prec, tree, f);
    }

    template <int D>
    void checkComplexTreeMatchesParts(const FunctionTree<D, double>& realTree,
                                      const FunctionTree<D, double>& imagTree,
                                      const FunctionTree<D, ComplexDouble>& ctree) {
        Coord<D> r1{};
        Coord<D> r2{};

        for (int i = 0; i < D; ++i) {
            r1[i] = 0.10 * (i + 1);
            r2[i] = -0.15 * (i + 1);
        }

        for (const auto& r : {r1, r2}) {
            auto z  = ctree.evalf(r);
            auto re = realTree.evalf(r);
            auto im = imagTree.evalf(r);

            REQUIRE(z.real() == Catch::Approx(re));
            REQUIRE(z.imag() == Catch::Approx(im));
        }
    }

    template <int D>
    void testComplexConstructor(bool useRealPart, bool useImagPart) {
        MultiResolutionAnalysis<D>* mra = nullptr;
        initialize(&mra);

        FunctionTree<D, double> realTree(*mra);
        FunctionTree<D, double> imagTree(*mra);

        Coord<D> posRe{};
        Coord<D> posIm{};
        std::array<int, static_cast<std::size_t>(D)> powerRe{};
        std::array<int, static_cast<std::size_t>(D)> powerIm{};

        for (int i = 0; i < D; ++i) {
            posRe[i] =  0.10 * (i + 1);
            posIm[i] = -0.07 * (i + 1);
            powerRe[i] = 0;
            powerIm[i] = (i == 0 ? 1 : 0);
        }

        if (useRealPart) {
            makeNonZeroTree<D>(realTree, 1.0, 0.4, posRe, powerRe);
        } else {
            realTree.setZero();
        }

        if (useImagPart) {
            makeNonZeroTree<D>(imagTree, 0.6, 0.7, posIm, powerIm);
        } else {
            imagTree.setZero();
        }

        FunctionTree<D, ComplexDouble> ctree(realTree, imagTree, nullptr, "complex_from_two_real");

        checkComplexTreeMatchesParts(realTree, imagTree, ctree);

        finalize(&mra);
    }


} // namespace function_tree_constructors