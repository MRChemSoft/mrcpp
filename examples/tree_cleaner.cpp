#include "MRCPP/MWFunctions"
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
    Printer::printHeader(0, "Testing TreeCleaner");

    // Constructing world box
    int corner[3] = {-1,-1,-1};
    int boxes[3]  = { 2, 2, 2};
    BoundingBox<3> world(min_scale, corner, boxes);

    // Constructing basis and MRA
    InterpolatingBasis basis(order);
    MultiResolutionAnalysis<3> MRA(world, basis);

    // Defining analytic function
    auto f = [] (const double *r) -> double {
        const double beta = 100.0;
        const double alpha = pow(beta/pi, 3.0/2.0);;
        double R = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        return alpha*exp(-beta*R*R);
    };

    // Setting up projector
    GridCleaner<3> clean(prec, max_scale);
    MWProjector<3> project(-1.0, max_scale);

    // Initializing MW function
    FunctionTree<3> f_tree(MRA);

    int iter = 0;
    int n_nodes = 1;
    while (n_nodes > 0) {
        project(f_tree, f); // Projecting on fixed grid
        n_nodes = clean(f_tree); // Refine grid and clear MW coefs 
        printout(0, " iter " << std::setw(3) << iter++ << std::setw(45));
        printout(0, " n_nodes " << std::setw(5) << n_nodes << std::endl);
    }
    // Projecting on final converged grid
    project(f_tree, f);

    double integral = f_tree.integrate();
    double sq_norm = f_tree.getSquareNorm();
    Printer::printSeparator(0, '-');
    Printer::printDouble(0, "Integral", integral);
    Printer::printDouble(0, "Square norm", sq_norm);

    timer.stop();
    Printer::printFooter(0, timer, 2);

    return 0;
}
