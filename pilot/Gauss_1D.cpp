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


    // Set parameters of the Gaussian
    double alpha = 0.5;
    double beta = 0.005;
    
    // Set center of the Gaussian
    double x0 = 0;
    double y0 = 0; // extra
    double z0 = 0; // extra

    // Create useful arrays 
    mrcpp::Coord<1> pos = {x0}; // -> this is the array holding the central position of the Gaussian
    std::array<int,1> pow = {0}; // -> this is the array holding the exponent of x prefactore of the Gaussian
    
    // Finally, we create the gaussian as a Gaussian function object, dim=1 with the parameters we set
    mrcpp::GaussFunc<1> my_gaussian(beta,alpha, pos, pow);

    // MRA - MW functions, first the BoundingBox that holds the limits of the box
    mrcpp::BoundingBox<1> bb({-5,5});

    // Then we create the MultiResolutionAnalysis object
    int order = 5;
    int MaxLevel = 10;

    mrcpp::MultiResolutionAnalysis mra(bb,order, MaxLevel);

   // Create the tree
   mrcpp::FunctionTree<1> f_tree(mra);

   //mrcpp::RepresentableFunction<1,double> my_gaussian_rep(my_gaussian, tree);
   mrcpp::project(0.01,f_tree, my_gaussian);



    // Plot the function
    mrcpp::Plotter<1> plot({-5});
    plot.setRange({10});
    plot.linePlot({1000}, f_tree, "F_tree");
    plot.linePlot({100}, my_gaussian, "Gaussian_1D");






    print::footer(0, timer, 2);
    return 0;
}