#include <cstdlib>
#include <numeric>
#include <iostream>
#include <vector>


#include "../api/Printer"
#include "../api/Timer"
#include "../api/Gaussians"
#include "../api/Plotter"
#include "../api/MWFunctions"
#include "../api/MWOperators"



using namespace mrcpp;

// DEBUG
bool debug = true;

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

#include "my_utilities.h"
#include "ZORA_utilities.h"

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
    mrcpp::CompFunction<3> VPhi_tree_top(mra);  // -> This is the tree that will hold the function
    mrcpp::CompFunction<3> VPhi_tree_bottom(mra);  // -> This is the tree that will hold the function
    mrcpp::CompFunction<3> GVPhi_tree_top(mra); // -> This is the tree that will hold the output of the convolution
    mrcpp::CompFunction<3> GVPhi_tree_bottom(mra); // -> This is the tree that will hold the output of the convolution
    mrcpp::CompFunction<3> core_el_tree(mra); // -> This is the tree that will hold the [nucleus-electron] potential

    mrcpp::CompFunction<3> Potential_tree(mra); // -> This is the tree that will hold the total potential
    // Then the 2 auxiliary trees
    mrcpp::CompFunction<3> Normalized_difference_tree(mra);
    mrcpp::CompFunction<3> Phi_n_copy_tree(mra);




    
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
    mrcpp::GaussFunc<3> Phi_trial_top(beta, alpha, pos, pow);
    mrcpp::GaussFunc<3> Phi_trial_bottom(beta, alpha, pos, pow);



    // 4) Let's DEFINE now the POTENTIAL for CORE ELECTRON function V_ee(r) = -1/r
    std::function<double(const Coord<3> &x)> V_ce = [] (const mrcpp::Coord<3> &r) -> double {
        auto R = std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        return -Z/R;
    };

    // 5) PROJECT the function ON the TREE
    mrcpp::project(GVPhi_tree_top, Phi_trial_top, building_precision); // I project the trial-function on the tree
    mrcpp::project(GVPhi_tree_bottom, Phi_trial_bottom, building_precision); // I project the trial-function on the tree

    double Psi_norm = std::sqrt(GVPhi_tree_top.getSquareNorm() + GVPhi_tree_bottom.getSquareNorm());
    GVPhi_tree_top.rescale(1.0/Psi_norm);
    GVPhi_tree_bottom.rescale(1.0/Psi_norm);

    // Create a vector of CompFunction<3> to hold the components
    std::vector<mrcpp::CompFunction<3>> Psi_2c;
    Psi_2c.push_back(GVPhi_tree_top);
    Psi_2c.push_back(GVPhi_tree_bottom);

    
    // As well as the potential:
    
    mrcpp::project(core_el_tree, V_ce, building_precision); 
    deep_copy(&Potential_tree, core_el_tree);
    
    if (n_electrons != 1){
        std::cout << "The number of elecrons is not supported" << '\n';
        return 0;
    }




    if (debug){
    // Befofe starting the SCF cycle, let's print some debugging information:
    std::cout << "Here is some debugging information: " << '\n';
    std::cout << "************************************" << '\n';
    std::cout << "Potential_tree Norm = ";
    std::cout << Potential_tree.getSquareNorm() << '\n';
    std::cout << "************************************" << '\n';
    std::cout << '\n'<< " Psi_trial (top) square norm = " << Psi_2c[0].getSquareNorm() << '\n';
    std::cout << '\n'<< " Psi_trial (bottom) square norm = " << Psi_2c[1].getSquareNorm() << '\n' << '\n';
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
    double apply_precision = building_precision; // Precision of the convolution

    //int maxIter = -1; //          -> Max number of iterations for the convolution, if -1 it goes untill the precision is reached
    
    double norm_diff = 1;   //    -> Norm of the difference in the 2 consecutive iterations (initialized as 1 to begin the while loop)

    int num_cycle =  0; //         -> SCF cycle counter 

    double E;


    // We now define all the trees that will be used to compute the enrgy and the SCF cycle
    mrcpp::CompFunction<3> K_tree(mra);  // -> This is the tree that will hold the K function
    mrcpp::CompFunction<3> K_inverted_tree(mra);  // -> This is the tree that will hold the K^-1 function
    std::vector<mrcpp::CompFunction<3> *> Nabla_K_tree(3); // -> This is the tree that will hold the gradient of K

    // K(r)
    std::function<double(const Coord<3> &x)> K_r = [] (const mrcpp::Coord<3> &r) -> double {
        double abs_r = std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        double constant = 1.0/(2.0 * m * c * c);
        return  abs_r / (abs_r - constant);
    };
    // K^-1(r)
    std::function<double(const Coord<3> &x)> K_inverted_r = [] (const mrcpp::Coord<3> &r) -> double {
        double abs_r = std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        double constant = 1.0/(2.0 * m * c * c);
        return (1 - (constant / abs_r));
    };

    // Project the K and K^-1 functions on the tree
    mrcpp::project(K_tree, K_r,building_precision);
    mrcpp::project(K_inverted_tree, K_inverted_r,building_precision);

    // Gradient of K
    mrcpp::ABGVOperator<3> D(mra, 0.0, 0.0); // deine the ABGV operator
    Nabla_K_tree = mrcpp::gradient(D,K_tree); // Gradient of K

    // Now we compute the Gradient of Psi_2c as well
    std::vector<mrcpp::CompFunction<3>*> Nabla_Psi_t;
    std::vector<mrcpp::CompFunction<3>*> Nabla_Psi_b;  
    Nabla_Psi_t = mrcpp::gradient(D, Psi_2c[0]);
    Nabla_Psi_b = mrcpp::gradient(D, Psi_2c[1]);

    // Create a vector to hold both Nabla_Psi_t and Nabla_Psi_b
    std::vector<std::vector<mrcpp::CompFunction<3>*>> Nabla_Psi_2c = {Nabla_Psi_t, Nabla_Psi_b};













/*     // Print the norm of the difference
    std::cout << "Cycle " << num_cycle << " done...  Norm of the difference = "<< norm_diff << '\n';
    E = energy_ZORA(mra, Psi_2c, Potential_tree);
    std::cout << "Energy = " << E << '\n';
    
    //std::cout << "Energy = " << E << '\n';
    mu = std::sqrt(-2.0 * E);
    std::cout << "mu = " << mu << '\n';
    std::cout << "************************************************************" << '\n';
    num_cycle++;
 */

    

    /* while (norm_diff > epsilon) {
        


        // Clear the trees
        VPhi_tree.clear();
        Phi_n_copy_tree.clear();
        Normalized_difference_tree.clear();

        // Phi_n_copy_tree = GVPhi_tree (If it's the first step then GVPhi_tree is the trial function)
        GVPhi_tree.deep_copy(&Phi_n_copy_tree); // deep copy copies the content of the first into the second

        // VPhi_tree = -2 * Potential_tree * GVPhi_tree
        mrcpp::multiply(apply_precision, VPhi_tree, -2.0, Potential_tree, Phi_n_copy_tree); 
        

        // Initialize GVPhi_tree
        GVPhi_tree.clear();

        // GVPhi_tree = \Tilde{\Phi^{n+1}}= - 2 * G^\mu (V \Phi^{n})
                                                                                                //mrcpp::apply(apply_precision, GVPhi_tree, Helm, VPhi_tree, maxIter, false);
        
        apply_Helmholtz(mu,building_precision, mra, GVPhi_tree, VPhi_tree);

        //std::cout << "Energy update: " << energy_update(mra, GVPhi_tree, Phi_n_copy_tree, Potential_tree) << '\n';
        

        // Normalize the function: \Tilde{\Phi^{n+1}} --> \Phi^{n+1}
        GVPhi_tree.normalize();
    
        // Normalized_difference_tree = \Phi^{n+1} - \Phi^{n}
        mrcpp::add(apply_precision, Normalized_difference_tree, +1.0, GVPhi_tree, -1.0, Phi_n_copy_tree);
        // norm_diff = || \Phi^{n+1} - \Phi^{n} ||
        norm_diff = Normalized_difference_tree.getSquareNorm();
        norm_diff = std::sqrt(norm_diff);





        // Print the norm of the difference
        std::cout << "Cycle " << num_cycle << " done...  Norm of the difference = "<< norm_diff << '\n';
        
        E = energy_ZORA(building_precision, mra, GVPhi_tree, Potential_tree);
        mu = std::sqrt(-2.0 * E);
        std::cout << "mu = " << mu << '\n';
        std::cout << '\n';
        std::cout << "************************************************************" << '\n';

        // Let's print all the steps!
        //plot.surfPlot({100,100}, GVPhi_tree, "Output_3D");
        //system("cat Output_3D.surf >> SCF_3D_Animation.surf");
        
        
        // Increment the cycle counter
        num_cycle++;
    }
 */
    print::footer(0, timer, 2);
    return 0;
}