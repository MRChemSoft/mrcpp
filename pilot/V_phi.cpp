#include "MRCPP/Printer"
#include "MRCPP/Timer"
#include "MRCPP/Gaussians"
#include "MRCPP/Plotter"
#include "MRCPP/MWFunctions"
// This is your personal sandbox to test out ideas.
// See the examples/ directory for inspiration.
// Please do not commit your pilot code to git.

using namespace mrcpp;

int main(int argc, char **argv) {
    Timer timer;

    // Initialize printing
    int printlevel = 0;
    Printer::init(printlevel);
    print::environment(0);
    print::header(0, "MRCPP pilot code");



    // We set the max order and max level of the MRA
    int order = 6;
    int MaxLevel = 15;

    std::array<int,2> box = {-5,5};
    mrcpp::BoundingBox<3> bb(box);
    mrcpp::MultiResolutionAnalysis mra(bb,order, MaxLevel);
    // Also, we create the tree
    mrcpp::FunctionTree<3> f_tree(mra);

    // My function will be now a 1s H atom orbital in 3D, given by the following lambda function




    std::function<double(const Coord<3> &x)> f = [](const Coord<3> &x) -> double {
        double r = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
        r = std::sqrt(r);
        double abs_r = std::abs(r);
        // In atomic units the Bohr radius is 1
        return std::exp(-abs_r);
    };

    // Let's define now the potential function V(r) = 1/r
    std::function<double(const Coord<3> &x)> V = [](const Coord<3> &x) -> double {
        double r = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
        r = std::sqrt(r);
        double abs_r = std::abs(r);
        // Prevent the singularity at r=0
        if(abs_r < 100) return -1.0 / abs_r;
        return -100;
        
    };

    // Now I define a whole new lambda function that will be the product of the two functions
    std::function<double(const Coord<3> &x)> f_V = [&](const Coord<3> &x) -> double {
        return f(x) * V(x);
    };



    // Project the function into the tree
    mrcpp::project<3>(0.0001,f_tree, f_V);

    // Plot the function
    mrcpp::Plotter<3> plot({-5,-5,0});
    plot.setRange({10,0,0},{0,10,0});
    plot.surfPlot({100,100}, f_tree, "V_phi_3D");
    //plot.linePlot({500}, f, "r-1");



    print::footer(0, timer, 2);
    return 0;
}