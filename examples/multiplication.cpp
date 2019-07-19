#include "MRCPP/Gaussians"
#include "MRCPP/MWFunctions"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include <array>

const auto min_scale = -4;
const auto max_depth = 25;

const auto order = 7;
const auto prec = 1.0e-5;

const auto D = 3;

int main(int argc, char **argv) {
    auto timer = mrcpp::Timer();

    // Initialize printing
    auto printlevel = 0;
    mrcpp::Printer::init(printlevel);
    mrcpp::print::environment(0);
    mrcpp::print::header(0, "Multiplying MW functions");

    // Constructing world box
    auto corner = std::array<int, D>{-1, -1, -1};
    auto boxes = std::array<int, D>{2, 2, 2};
    auto world = mrcpp::BoundingBox<D>(min_scale, corner, boxes);

    // Constructing basis and MRA
    auto basis = mrcpp::InterpolatingBasis(order);
    auto MRA = mrcpp::MultiResolutionAnalysis<D>(world, basis, max_depth);

    // Setting up analytic Gaussians
    auto beta = 20.0;
    auto alpha = std::pow(beta / mrcpp::pi, 3.0 / 2.0);
    auto f_pos = mrcpp::Coord<D>{0.0, 0.0, 0.1};
    auto g_pos = mrcpp::Coord<D>{0.0, 0.0, -0.1};
    auto power = std::array<int, D>{0, 0, 0};
    auto f_func = mrcpp::GaussFunc<D>(beta, alpha, f_pos, power);
    auto g_func = mrcpp::GaussFunc<D>(beta, alpha, g_pos, power);

    // Initialize MW functions
    mrcpp::FunctionTree<D> f_tree(MRA);
    mrcpp::FunctionTree<D> g_tree(MRA);
    mrcpp::FunctionTree<D> h_tree(MRA);

    // Projecting f and g
    mrcpp::project<D>(prec, f_tree, f_func);
    mrcpp::project<D>(prec, g_tree, g_func);

    // h = f*g
    mrcpp::multiply(prec, h_tree, 1.0, f_tree, g_tree);

    auto integral = h_tree.integrate();
    auto sq_norm = h_tree.getSquareNorm();
    mrcpp::print::value(0, "Integral", integral);
    mrcpp::print::value(0, "Square norm", sq_norm);
    mrcpp::print::footer(0, timer, 2);

    return 0;
}
