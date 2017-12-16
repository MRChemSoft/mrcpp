#include "MRCPP/MWFunctions"
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
    Printer::printHeader(0, "Multiplying MW functions");

    // Constructing world box
    int corner[3] = {-1,-1,-1};
    int boxes[3]  = { 2, 2, 2};
    BoundingBox<3> world(min_scale, corner, boxes);

    // Constructing basis and MRA
    InterpolatingBasis basis(order);
    MultiResolutionAnalysis<3> MRA(world, basis);

    // Setting up projectors
    MWProjector<3> project(prec, max_scale);
    MWMultiplier<3> mult(prec, max_scale);

    // Setting up analytic Gaussians
    double beta = 20.0;
    double alpha = pow(beta/pi, 3.0/2.0);
    double f_pos[3] = {0.0, 0.0,  0.1};
    double g_pos[3] = {0.0, 0.0, -0.1};
    GaussFunc<3> f_func(beta, alpha, f_pos);
    GaussFunc<3> g_func(beta, alpha, g_pos);

    // Initialize MW functions
    FunctionTree<3> f_tree(MRA);
    FunctionTree<3> g_tree(MRA);
    FunctionTree<3> h_tree(MRA);

    // Projecting f and g
    project(f_tree, f_func);
    project(g_tree, g_func);

    // h = f*g
    mult(h_tree, 1.0, f_tree, g_tree);

    double integral = h_tree.integrate();
    double sq_norm = h_tree.getSquareNorm();
    Printer::printDouble(0, "Integral", integral);
    Printer::printDouble(0, "Square norm", sq_norm);

    timer.stop();
    Printer::printFooter(0, timer, 2);

    return 0;
}

