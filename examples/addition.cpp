#include <tuple>

#include "MRCPP/Gaussians"
#include "MRCPP/MWFunctions"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

const auto min_scale = -4;
const auto max_depth = 25;

const auto order = 7;
const auto prec = 1.0e-5;
const auto D = 3; // Dimensions

int main(int argc, char **argv) {
    auto timer = mrcpp::Timer();

    // Initialize printing
    auto printlevel = 0;
    mrcpp::Printer::init(printlevel);
    mrcpp::print::environment(0);
    mrcpp::print::header(0, "Adding MW functions");

    // Constructing world box
    auto corner = std::array<int, D>{-1, -1, -1};
    auto boxes = std::array<int, D>{2, 2, 2};
    auto world = mrcpp::BoundingBox<D>(min_scale, corner, boxes);

    // Constructing basis and MRA
    auto basis = mrcpp::InterpolatingBasis(order);
    auto MRA = mrcpp::MultiResolutionAnalysis<3>(world, basis, max_depth);

    // Setting up analytic Gaussians
    auto beta = 20.0;
    auto alpha = std::pow(beta / mrcpp::pi, 3.0 / 2.0);
    auto pos_1 = std::array<double, D>{0.0, 0.0, 0.1};
    auto pos_2 = std::array<double, D>{0.0, 0.0, -0.1};
    auto pos_3 = std::array<double, D>{0.0, 0.0, 0.3};
    auto power = std::array<int, D>{0, 0, 0};
    auto f_func_1 = mrcpp::GaussFunc<D>(beta, alpha, pos_1, power);
    auto f_func_2 = mrcpp::GaussFunc<D>(beta, alpha, pos_2, power);
    auto f_func_3 = mrcpp::GaussFunc<D>(beta, alpha, pos_3, power);

    // Initialize MW functions
    mrcpp::FunctionTree<D> f_tree_1(MRA);
    mrcpp::FunctionTree<D> f_tree_2(MRA);
    mrcpp::FunctionTree<D> f_tree_3(MRA);
    mrcpp::FunctionTree<D> g_tree(MRA);

    // Projecting functions
    mrcpp::project(prec, f_tree_1, f_func_1);
    mrcpp::project(prec, f_tree_2, f_func_2);
    mrcpp::project(prec, f_tree_3, f_func_3);

    // g = f_1 - 2*f_2 + 3*f_3
    auto f_vec = mrcpp::FunctionTreeVector<D>();
    f_vec.push_back(std::make_tuple(1.0, &f_tree_1));
    f_vec.push_back(std::make_tuple(-2.0, &f_tree_2));
    f_vec.push_back(std::make_tuple(3.0, &f_tree_3));
    mrcpp::add(prec, g_tree, f_vec);

    auto integral = g_tree.integrate();
    auto sq_norm = g_tree.getSquareNorm();
    mrcpp::print::value(0, "Integral", integral);
    mrcpp::print::value(0, "Square norm", sq_norm);
    mrcpp::print::footer(0, timer, 2);

    return 0;
}
