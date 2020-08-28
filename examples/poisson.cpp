#include "MRCPP/Gaussians"
#include "MRCPP/MWFunctions"
#include "MRCPP/MWOperators"
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
    mrcpp::print::environment(0);
    mrcpp::print::header(0, "Applying Poisson operator");

    // Constructing world box
    auto corner = std::array<int, D>{-1, -1, -1};
    auto boxes = std::array<int, D>{2, 2, 2};
    auto world = mrcpp::BoundingBox<D>(min_scale, corner, boxes);

    // Constructing basis and MRA
    auto basis = mrcpp::InterpolatingBasis(order);
    auto MRA = mrcpp::MultiResolutionAnalysis<D>(world, basis, max_depth);

    // Setting up analytic Gaussian
    auto beta = 100.0;
    auto alpha = std::pow(beta / mrcpp::pi, 3.0 / 2.0);
    auto pos = mrcpp::Coord<D>{mrcpp::pi / 3.0, mrcpp::pi / 3.0, mrcpp::pi / 3.0};
    auto power = std::array<int, D>{0, 0, 0};
    auto f_func = mrcpp::GaussFunc<D>(beta, alpha, pos, power);

    // Computing analytic energy
    auto ana_energy = f_func.calcCoulombEnergy(f_func);

    // Initializing MW functions and operator
    mrcpp::print::separator(0, ' ');
    mrcpp::print::memory(0, "used memory pre Poisson");
    mrcpp::PoissonOperator P(MRA, prec);
    mrcpp::print::memory(0, "used memory post Poisson");

    // Projecting function
    auto t1 = mrcpp::Timer();
    mrcpp::print::separator(0, ' ');
    mrcpp::print::memory(0, "used memory pre project");
    mrcpp::FunctionTree<D> f_tree(MRA);
    mrcpp::build_grid(f_tree, f_func);
    mrcpp::project(prec, f_tree, f_func);
    mrcpp::print::memory(0, "used memory post project");
    t1.stop();

    // Applying Poisson operator
    auto t2 = mrcpp::Timer();
    mrcpp::print::separator(0, ' ');
    mrcpp::print::memory(0, "used memory pre apply");
    mrcpp::FunctionTree<D> g_tree(MRA);
    mrcpp::apply(prec, g_tree, P, f_tree);
    mrcpp::print::memory(0, "used memory post apply");
    t2.stop();

    auto f_int = f_tree.integrate();
    auto f_norm = std::sqrt(f_tree.getSquareNorm());
    auto g_int = g_tree.integrate();
    auto g_norm = std::sqrt(g_tree.getSquareNorm());
    auto num_energy = mrcpp::dot(g_tree, f_tree);
    auto error = (num_energy - ana_energy) / num_energy;

    mrcpp::print::separator(0, ' ');
    mrcpp::print::tree(0, "f_tree", f_tree, t1);
    mrcpp::print::tree(0, "g_tree", g_tree, t2);
    mrcpp::print::separator(0, ' ');
    mrcpp::print::value(0, "f_tree integral", f_int);
    mrcpp::print::value(0, "f_tree norm", f_norm);
    mrcpp::print::separator(0, ' ');
    mrcpp::print::value(0, "g_tree integral", g_int);
    mrcpp::print::value(0, "g_tree norm", g_norm);
    mrcpp::print::separator(0, ' ');
    mrcpp::print::value(0, "Analytic energy", ana_energy, "(au)");
    mrcpp::print::value(0, "Numerical energy", num_energy, "(au)");
    mrcpp::print::value(0, "Error", error, "(rel)", 1);
    mrcpp::print::separator(0, ' ');
    mrcpp::print::footer(0, timer, 2);

    return 0;
}
