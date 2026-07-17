#include "MRCPP/Printer"
#include "MRCPP/Timer"
#include "MRCPP/Gaussians"
#include "MRCPP/Plotter"
#include "MRCPP/MWFunctions"
#include "MRCPP/MWOperators"

#include <cstdlib>
#include <numeric>

using namespace mrcpp;

// DEBUG
bool debug = false;

// GLOBAL VARIABLES
int Z;
double c = 137.035999139;
double m = 1.0;
int n_electrons;


int order;
int MaxLevel;
double building_precision;
double epsilon;
int Relativity;

/*
 * ==================================================================================================================================
 * This program will introduce relativistic effects in the H atom, in detail, with the ZORA approximation.
 * 
 * In this framework we can say that the Hamiltonian is given by:
 *      H_{ZORA} = T_{ZORA} + V
 * V Is the same as the H atom in general. While T_{ZORA} is different. From classical to relativistic:
 *      T = - p^2/2m -> T_{ZORA} = p \cdot K(r) p
 * Where this K term is given by:
 *     K(r) = [1-V/(2mc^2)]^{-1}
 * Remember that each p is an operator, so the product is not trivial.
 * ==================================================================================================================================
*/






#include "He_atom_utils.h"







// ==================================================================================================================================


// ==================================================================================================================================







// ==================================================================================================================================

// Funzione per leggere i parametri dal file
std::unordered_map<std::string, double> readParameters(const std::string &filename) {
    std::unordered_map<std::string, double> parameters;
    std::ifstream infile(filename);
    std::string line;

    while (std::getline(infile, line)) {
        std::istringstream iss(line);
        std::string key;
        double value;
        if (iss >> key >> value) {
            parameters[key] = value;
        }
    }

    return parameters;
}







// ==================================================================================================================================
// ================================================ START OF THE MAIN FUNCTION ======================================================
// ==================================================================================================================================

int main(int argc, char **argv) {
    Timer timer;

    // Initialize printing
    int printlevel = 3;
    Printer::init(printlevel);
    print::environment(0);
    print::header(0, "MRCPP pilot code");
    
// ==================================================================================================================================
// ================= Start of the parameters reading ================================================================================
// ==================================================================================================================================
    auto parameters = readParameters("Input_Parameters.txt");

    // Assegna i valori ai parametri corrispondenti
    order = static_cast<int>(parameters["order"]);
    MaxLevel = static_cast<int>(parameters["MaxLevel"]);
    building_precision = parameters["building_precision"];
    epsilon = parameters["epsilon"];
    Z = static_cast<int>(parameters["Z"]);
    // Set the charge as neutra by:
    Z = 2;
    n_electrons = Z;
    Relativity = static_cast<int>(parameters["Relativity"]);

    std::cout << "Parameters read from file: " << '\n' << '\n';
    std::cout << " Order =" << '\t'<< '\t' << order << '\n';
    std::cout << " MaxLevel =" << '\t'<< '\t' << MaxLevel << '\n';
    std::cout << " Building precision =" << '\t' << building_precision << '\n';
    std::cout << " Epsilon =" << '\t'<< '\t' << epsilon << '\n';
    std::cout << " Z =" << '\t'<< '\t' << '\t' << Z << '\n';
    std::cout << " Number of electrons =" << '\t' << n_electrons << '\n';
    if (Relativity == 0){
        std::cout << " Relativity ="<< '\t'<< '\t' <<"Non-Relativistic" << '\n';
    }
    else if (Relativity == 1){
        std::cout << " Relativity =" << '\t'<< '\t' << "ZORA" << '\n';
    }
    else{
        std::cout << " Relativity =" << '\t'<< '\t' << "not supported" << '\n';
        return 0;
    }
    std::cout << '\n' << '\n' << "************************************************************" << '\n';
    




    // ==================================================================================================================================
    // ================= End of the parameters reading ==================================================================================
    // ==================================================================================================================================

    // 1) Build the Box and MRA
    std::array<int,2> box = {-20,20};
    mrcpp::BoundingBox<3> bb(box);
    mrcpp::MultiResolutionAnalysis mra(bb, order, MaxLevel);
    // Also, we create the trees for each function we'll be dealing with
    mrcpp::FunctionTree<3> VPhi_tree(mra);  // -> This is the tree that will hold the function
    mrcpp::FunctionTree<3> Phi_in(mra); // -> This is the tree that will hold the output of the convolution
    mrcpp::FunctionTree<3> J_tree(mra); // -> This is the tree that will hold the [electron-electron] (pair) potential
    mrcpp::FunctionTree<3> Potential_tree(mra); // -> This is the tree that will hold the total potential
    // Then the 2 auxiliary trees
    mrcpp::FunctionTree<3> Normalized_difference_tree(mra);
    mrcpp::FunctionTree<3> Phi_n_copy_tree(mra);

    // ==================================================================================================================================
    // ================= My function will be now a 1s H atom orbital in 3D, given by the following lambda function ======================
    // ==================================================================================================================================

    // 2) We define declare operator G^\mu as the Helmltz Convolutional Operator:      [G^\mu]
    double mu = 1.0;



    // 3) Define a trial function as a Gaussian 
    //   The function is defined as \Phi(r) = \alpha\exp(-\beta r^2)
    // Parameters:    
    double beta = 0.5;
    double alpha = (beta/3.14159265358979)*(beta/3.14159265358979)*(beta/3.14159265358979);
    alpha = std::sqrt(alpha);
    mrcpp::Coord<3> pos = {0,0,0};
    std::array<int,3> pow = {0,0,0};
    // Gaussian function
    //mrcpp::GaussFunc<3> Phi_trial(beta, alpha, pos, pow);

    // Define the trial function as a Slater-type orbital
    std::function<double(const Coord<3> &x)> slater = [] (const mrcpp::Coord<3> &r) -> double {
        auto R = std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        double lambda = std::sqrt(1. - (1/(c*c)) );
        double alpha = 1.69;
        return exp(-alpha*R);
    };
    
    // Define another function that is zero everywhere
    std::function<double(const Coord<3> &x)> zero = [] (const mrcpp::Coord<3> &r) -> double {
        return 0.0;
    };
        

    // 4) Let's DEFINE now the POTENTIAL for CORE ELECTRON function V_ee(r) = -1/r
    std::function<double(const Coord<3> &x)> V_ce = [] (const mrcpp::Coord<3> &r) -> double {
        auto R = std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        return -Z/R;
    };






    // 5) PROJECT the function ON the TREE
    mrcpp::project<3>(building_precision, Phi_in, slater); // I project the trial-function on the tree
    mrcpp::project<3>(building_precision, Potential_tree, V_ce); // I project the potential on the tree

    Phi_in.normalize();
    //std::cout << "Psi_trial" << Phi_in << '\n';




    if (debug){
    // Befofe starting the SCF cycle, let's print some debugging information:
    std::cout << "Here is some debugging information: " << '\n';
    std::cout << "************************************" << '\n';
    std::cout << "Potential_tree = " << '\n';
    std::cout << Potential_tree << '\n';
    std::cout << "************************************" << '\n';
    std::cout << "Gauss_tree = " << '\n';
    std::cout << Phi_in << '\n';
    std::cout << "************************************" << '\n';
    }
    /*
     *-----------------------------------
     *-----------------------------------
     *         BEGIN SCF CYCLE
     *-----------------------------------
     *-----------------------------------
    */

    // Parameters for the SCF
    double norm_diff = 1;   //    -> Norm of the difference in the 2 consecutive iterations (initialized as 1 to begin the while loop)
    int num_cycle =  1; //         -> SCF cycle counter 
    double E;



    FunctionTree<3> G_V_Phi(mra);

 
    while (norm_diff > epsilon) {
        // Clear the trees
        std::cout << '\n';
        std::cout << "***********************************************" << '\n';

        Normalized_difference_tree.clear();

        compute_He_energy(mra, Phi_in, Potential_tree, mu, E);
        std::cout << "> mu = " << mu << '\n';

        G_V_Phi.clear();
        // Apply the Helmholtz operator to the trial function
        apply_Helmholtz_He(mu,mra, G_V_Phi, Phi_in, Potential_tree);

        // Normalize the result
        G_V_Phi.normalize();

        // Find the difference:
        mrcpp::add(building_precision, Normalized_difference_tree, 1.0, G_V_Phi, -1.0, Phi_in);

        norm_diff = Normalized_difference_tree.getSquareNorm();
        norm_diff = std::sqrt(norm_diff);



        // Deep copy the new function to the input function for the next iteration
        Phi_in.clear();
        G_V_Phi.deep_copy(&Phi_in);

        // Print the norm of the difference
        std::cout << "Cycle " << num_cycle << " done. --->  * Norm of the difference = "<< norm_diff << '\n';
        
        
        // Increment the cycle counter
        num_cycle++;
    }

    print::footer(0, timer, 2);
    return 0;
}