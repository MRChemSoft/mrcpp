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
    mrcpp::BoundingBox<1> bb({-5,5});
    mrcpp::MultiResolutionAnalysis mra(bb,order, MaxLevel);
    // Also, we create the tree
    mrcpp::FunctionTree<1> f_tree(mra);


    /*
        * Simples approach, we use the build in function of Gaussians
        * PROBLEM: The function is not normalized AND the function in not in absoulte value
    
    Define the Function as 1/r:
    
    mrcpp::Coord<1> pos = {0}; // -> we are working in 1D, so we only need one coordinate. We set the center at 0
    std::array<int,1> pow = {-1}; // In this way we define the function as 1/r
     * Instead of creatin a brand new class, 
     * we use the gaussian class, where the prefactor 
     * is 1 and the exponent (beta) is 0. 
     * Choosing the exponent of x as -1
    
    double alpha = 1;
    double beta = 0;
    mrcpp::GaussFunc<1> my_gaussian(beta,alpha, pos, pow);
    */    

    std::function<double(const Coord<1> &r)> f = [](const Coord<1> &r) -> double {
        double abs_r = std::abs(r[0]);
        double max  = 50;
        if(abs_r >= max){
        return -max;
        }
        else{
        return -1.0 / abs_r;
        }
    };






    // Project the function into the tree
    mrcpp::project<1>(0.0001,f_tree, f);

    // Plot the function
    mrcpp::Plotter<1> plot({-5});
    plot.setRange({10});
    plot.linePlot({1000}, f_tree, "P_r-1");
    //plot.linePlot({500}, f, "r-1");



    print::footer(0, timer, 2);
    return 0;
}