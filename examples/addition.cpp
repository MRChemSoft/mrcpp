#include "MRCPP/MWFunctions"
#include "MRCPP/Gaussians"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

const int min_scale = -4;
const int max_depth = 25;

const int order = 7;
const double prec = 1.0e-5;

int main(int argc, char **argv) {
    Timer timer;

    // Initialize printing
    int printlevel = 0;
    Printer::init(printlevel);
    Printer::printEnvironment();
    Printer::printHeader(0, "Adding MW functions");

    // Constructing world box
    int corner[3] = {-1,-1,-1};
    int boxes[3]  = { 2, 2, 2};
    BoundingBox<3> world(min_scale, corner, boxes);

    // Constructing basis and MRA
    InterpolatingBasis basis(order);
    MultiResolutionAnalysis<3> MRA(world, basis, max_depth);

    // Setting up analytic Gaussians
    double beta = 20.0;
    double alpha = pow(beta/pi, 3.0/2.0);
    double pos_1[3] = {0.0, 0.0,  0.1};
    double pos_2[3] = {0.0, 0.0, -0.1};
    double pos_3[3] = {0.0, 0.0,  0.3};
    GaussFunc<3> f_func_1(beta, alpha, pos_1);
    GaussFunc<3> f_func_2(beta, alpha, pos_2);
    GaussFunc<3> f_func_3(beta, alpha, pos_3);

    // Initialize MW functions
    FunctionTree<3> f_tree_1(MRA);
    FunctionTree<3> f_tree_2(MRA);
    FunctionTree<3> f_tree_3(MRA);
    FunctionTree<3> g_tree(MRA);

    // Projecting functions
    mrcpp::project(prec, f_tree_1, f_func_1);
    mrcpp::project(prec, f_tree_2, f_func_2);
    mrcpp::project(prec, f_tree_3, f_func_3);

    // g = f_1 - 2*f_2 + 3*f_3
    FunctionTreeVector<3> f_vec;
    f_vec.push_back( 1.0, &f_tree_1);
    f_vec.push_back(-2.0, &f_tree_2);
    f_vec.push_back( 3.0, &f_tree_3);
    mrcpp::add(prec, g_tree, f_vec);

    double integral = g_tree.integrate();
    double sq_norm = g_tree.getSquareNorm();
    Printer::printDouble(0, "Integral", integral);
    Printer::printDouble(0, "Square norm", sq_norm);

    timer.stop();
    Printer::printFooter(0, timer, 2);

    return 0;
}

