#include "MRCPP/MWFunctions"
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
    mrcpp::Printer::printHeader(0, "Projecting analytic function");

    // Constructing world box
    int corner[3] = {-1,-1,-1};
    int boxes[3]  = { 2, 2, 2};
    mrcpp::BoundingBox<3> world(min_scale, corner, boxes);

    // Constructing basis and MRA
    mrcpp::InterpolatingBasis basis(order);
    mrcpp::MultiResolutionAnalysis<3> MRA(world, basis, max_depth);

    // Defining analytic function
    auto f = [] (const double *r) -> double {
        const double beta = 100.0;
        const double alpha = pow(beta/mrcpp::pi, 3.0/2.0);;
        double R = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        return alpha*exp(-beta*R*R);
    };

    // Projecting function
    mrcpp::FunctionTree<3> f_tree(MRA);
    mrcpp::project(prec, f_tree, f);

    double integral = f_tree.integrate();
    double sq_norm = f_tree.getSquareNorm();
    mrcpp::Printer::printDouble(0, "Integral", integral);
    mrcpp::Printer::printDouble(0, "Square norm", sq_norm);

    timer.stop();
    mrcpp::Printer::printFooter(0, timer, 2);

    return 0;
}

