#include "MRCPP/MWFunctions"
#include "MRCPP/MWOperators"
#include "MRCPP/Gaussians"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

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
    mrcpp::Printer::printEnvironment();
    mrcpp::Printer::printHeader(0, "Applying Poisson operator");

    // Constructing world box
    auto corner = std::array<int, D>{-1,-1,-1};
    auto boxes  = std::array<int, D>{2, 2, 2};
    auto world = mrcpp::BoundingBox<D>(min_scale, corner, boxes);

    // Constructing basis and MRA
    auto basis = mrcpp::InterpolatingBasis(order);
    auto MRA = mrcpp::MultiResolutionAnalysis<D>(world, basis, max_depth);

    // Setting up analytic Gaussian
    auto beta = 100.0;
    auto alpha = std::pow(beta/mrcpp::pi, 3.0/2.0);
    auto pos = mrcpp::Coord<D>{mrcpp::pi/3.0, mrcpp::pi/3.0, mrcpp::pi/3.0};
    auto power = std::array<int, D>{0, 0, 0};
    auto f_func = mrcpp::GaussFunc<D>(beta, alpha, pos, power);

    // Computing analytic energy
    auto ana_energy = f_func.calcCoulombEnergy(f_func);

    // Initializing MW functions and operator
    auto P = mrcpp::PoissonOperator(MRA, prec);
    auto f_tree = mrcpp::FunctionTree<D>(MRA);
    auto g_tree = mrcpp::FunctionTree<D>(MRA);

    // Projecting function
    mrcpp::build_grid(f_tree, f_func);
    mrcpp::project(prec, f_tree, f_func);

    // Applying Poisson operator
    mrcpp::apply(prec, g_tree, P, f_tree);

    auto f_int = f_tree.integrate();
    auto f_norm = std::sqrt(f_tree.getSquareNorm());
    auto g_int = g_tree.integrate();
    auto g_norm = std::sqrt(g_tree.getSquareNorm());
    auto num_energy = mrcpp::dot(g_tree, f_tree);
    auto error = (num_energy-ana_energy)/num_energy;

    mrcpp::Printer::printSeparator(0, ' ');
    mrcpp::Printer::printDouble(0, "f_tree integral", f_int);
    mrcpp::Printer::printDouble(0, "f_tree norm", f_norm);
    mrcpp::Printer::printSeparator(0, ' ');
    mrcpp::Printer::printDouble(0, "g_tree integral", g_int);
    mrcpp::Printer::printDouble(0, "g_tree norm", g_norm);
    mrcpp::Printer::printSeparator(0, ' ');
    mrcpp::Printer::printDouble(0, "Analytic energy", ana_energy);
    mrcpp::Printer::printDouble(0, "Numerical energy", num_energy);
    mrcpp::Printer::printDouble(0, "Relative error", error, 1);
    mrcpp::Printer::printSeparator(0, ' ');

    timer.stop();
    mrcpp::Printer::printFooter(0, timer, 2);

    return 0;
}
