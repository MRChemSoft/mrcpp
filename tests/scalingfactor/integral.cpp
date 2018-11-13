#include "catch.hpp"

#include "factory_functions.h"

#include "MRCPP/MWFunctions"
#include "MRCPP/Gaussians"

using namespace mrcpp;

namespace scaling {

template<int D> void testScaling();

template<int D>
double alpha_gen(double sigma) {
    return std::pow(1.0/std::sqrt(2.0*pi*std::pow(sigma, 2.0)), D);
}

double beta_gen(double scaling_factor, double sigma) {
    return std::pow(scaling_factor, 2.0)/(2.0*std::pow(sigma, 2.0));
}

template<int D>
std::array<int, D> generate_corner() {
    auto corner = std::array<int, D> {};
    for (auto & x : corner)
        x = 0;
    return corner;
}

template<int D>
std::array<int, D> generate_boxes() {
    auto boxes = std::array<int, D> {};
    for (auto & x : boxes)
        x = 1;
    return boxes;
}

template<int D>
std::array<double, D> generate_scaling_factor(double sf) {
    auto scaling_factor = std::array<double, D> {};
    for (auto & x : scaling_factor)
        x = sf;
    return scaling_factor;
}

template<int D>
std::array<double, D> generate_pos(double x0) {
    auto pos = std::array<double, D> {};
    for (auto & x : pos) {
        x = x0;
    }
    return pos;
}

SCENARIO("Testing Function Values", "[scaling], [integral]") {
    GIVEN("1D") {
        testScaling<1>();
    }
    GIVEN("2D") {
     testScaling<2>();
    }
    GIVEN("3D") {
     testScaling<3>();
    }
}

template<int D> void testScaling() {
    const auto prec = 1.0e-3;
    auto min_scale = 0;

    auto corner = generate_corner<D>();
    auto boxes = generate_boxes<D>();
    auto sf = generate_scaling_factor<D>(2.0*pi);

    auto world = BoundingBox<D>(min_scale, corner, boxes, sf);
    auto basis = InterpolatingBasis(5);
    auto MRA = MultiResolutionAnalysis<D>(world, basis, 25);

    auto sigma = 0.2;
    auto alpha = alpha_gen<D>(sigma);
    auto beta = beta_gen(sf[0], sigma);
    auto power = std::array<int, D>{};
    auto pos = generate_pos<D>(0.5);
    auto gauss = GaussFunc<D>(beta, alpha, pos, power);

    FunctionTree<D> f_tree(MRA);
    project(prec, f_tree, gauss);

    WHEN("Scalings") {
        THEN("test") {
            REQUIRE(Approx(1.0) == f_tree.integrate() );
        }
    }

}

} // namespace scaling
