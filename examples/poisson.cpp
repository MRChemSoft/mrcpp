#include "MRCPP/MWFunctions"
#include "MRCPP/MWOperators"
#include "MRCPP/Gaussians"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

const int min_scale = -4;
const int max_depth = 25;

const int order = 7;
const double prec = 1.0e-5;

int main(int argc, char **argv) {
    mrcpp::Timer timer;

    // Initialize printing
    int printlevel = 0;
    mrcpp::Printer::init(printlevel);
    mrcpp::Printer::printEnvironment();
    mrcpp::Printer::printHeader(0, "Applying Poisson operator");

    // Constructing world box
    int corner[3] = {-1,-1,-1};
    int boxes[3]  = { 2, 2, 2};
    mrcpp::BoundingBox<3> world(min_scale, corner, boxes);

    // Constructing basis and MRA
    mrcpp::InterpolatingBasis basis(order);
    mrcpp::MultiResolutionAnalysis<3> MRA(world, basis, max_depth);

    // Setting up analytic Gaussian
    double beta = 100.0;
    double alpha = pow(beta/mrcpp::pi, 3.0/2.0);
    double pos[3] = {mrcpp::pi/3.0,mrcpp::pi/3.0,mrcpp::pi/3.0};
    mrcpp::GaussFunc<3> f_func(beta, alpha, pos);

    // Computing analytic energy
    double ana_energy = f_func.calcCoulombEnergy(f_func);

    // Initializing MW functions and operator
    mrcpp::PoissonOperator P(MRA, prec);
    mrcpp::FunctionTree<3> f_tree(MRA);
    mrcpp::FunctionTree<3> g_tree(MRA);

    // Projecting function
    mrcpp::build_grid(f_tree, f_func);
    mrcpp::project(prec, f_tree, f_func);

    // Applying Poisson operator
    mrcpp::apply(prec, g_tree, P, f_tree);

    double f_int = f_tree.integrate();
    double f_norm = sqrt(f_tree.getSquareNorm());
    double g_int = g_tree.integrate();
    double g_norm = sqrt(g_tree.getSquareNorm());
    double num_energy = g_tree.dot(f_tree);
    double error = (num_energy-ana_energy)/num_energy;

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

