#include "MRCPP/MWFunctions"
#include "MRCPP/MWOperators"
#include "MRCPP/Gaussians"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

const int min_scale = -4;
const int max_scale = 20;

const int order = 7;
const double prec = 1.0e-5;

int main(int argc, char **argv) {
    Timer timer;

    // Initialize printing
    int printlevel = 0;
    Printer::init(printlevel);
    Printer::printEnvironment();
    Printer::printHeader(0, "Applying Poisson operator");

    // Constructing world box
    int corner[3] = {-1,-1,-1};
    int boxes[3]  = { 2, 2, 2};
    BoundingBox<3> world(min_scale, corner, boxes);

    // Constructing basis and MRA
    InterpolatingBasis basis(order);
    MultiResolutionAnalysis<3> MRA(world, basis);

    // Setting up projectors
    GridGenerator<3> grid(max_scale);
    MWProjector<3> project(prec, max_scale);
    MWConvolution<3> apply(prec, max_scale);

    // Setting up analytic Gaussian
    double beta = 100.0;
    double alpha = pow(beta/pi, 3.0/2.0);
    double pos[3] = {pi/3.0,pi/3.0,pi/3.0};
    GaussFunc<3> f_func(beta, alpha, pos);

    // Computing analytic energy
    double ana_energy = f_func.calcCoulombEnergy(f_func);

    // Initializing MW functions and operator
    PoissonOperator P(MRA, prec);
    FunctionTree<3> f_tree(MRA);
    FunctionTree<3> g_tree(MRA);

    // Projecting function
    grid(f_tree, f_func);
    project(f_tree, f_func);

    // Applying Poisson operator
    apply(g_tree, P, f_tree);

    double f_int = f_tree.integrate();
    double f_norm = sqrt(f_tree.getSquareNorm());
    double g_int = g_tree.integrate();
    double g_norm = sqrt(g_tree.getSquareNorm());
    double num_energy = g_tree.dot(f_tree);
    double error = (num_energy-ana_energy)/num_energy;

    Printer::printSeparator(0, ' ');
    Printer::printDouble(0, "f_tree integral", f_int);
    Printer::printDouble(0, "f_tree norm", f_norm);
    Printer::printSeparator(0, ' ');
    Printer::printDouble(0, "g_tree integral", g_int);
    Printer::printDouble(0, "g_tree norm", g_norm);
    Printer::printSeparator(0, ' ');
    Printer::printDouble(0, "Analytic energy", ana_energy);
    Printer::printDouble(0, "Numerical energy", num_energy);
    Printer::printDouble(0, "Relative error", error, 1);
    Printer::printSeparator(0, ' ');

    timer.stop();
    Printer::printFooter(0, timer, 2);

    return 0;
}

