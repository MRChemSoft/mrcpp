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
#include <fstream>
#include <iostream>

void print_memory_usage() {
    std::ifstream statm("/proc/self/statm");
    long size, resident;
    statm >> size >> resident;
    std::cout << "Virtual Memory: " << 0.000001*size * 4 << " GB\t";
    std::cout << "Resident Set Size: " << 0.000001*resident * 4 << " GB\n";
}

using namespace mrcpp;
using Orbital_Pointer_List = std::vector<std::vector<mrcpp::CompFunction<3>>*>;
using Spinor_CompFunction = std::vector<mrcpp::CompFunction<3>>;
using Spinor_gradients_pointer_list = std::vector<std::vector<std::vector<mrcpp::CompFunction<3>*>>>;







// DEBUG
bool debug = false;
bool verbose = false;

// GLOBAL VARIABLES
int Z =4;
//double c = 137.035999146; // --> This is the speed of light in atomic units
double c = 137.035999177; 
double m = 1.0;
int n_electrons;


int order;
int MaxLevel;
double building_precision;
double epsilon;
int Relativity;
int num_cycle =  0;     //    -> SCF cycle counter 

// Define another function that is zero everywhere
static std::function<double(const Coord<3> &x)> zero = [] (const mrcpp::Coord<3> &r) -> double {
    return 0.0;
};
 



#include "my_utilities.h"
//#include "ZORA_utilities.h"

#include "Be_atom_ZORA_utils.h"
#include "Be_Update_var.h"
#include "Smeared_potential.h"
//#include "ZORA_utilities.h"







void SCF_Cycle_ZORA(MultiResolutionAnalysis<3> &MRA, 
    Orbital_Pointer_List &Spin_orbitals, 
    mrcpp::CompFunction<3> &core_el_tree, 
    std::vector<std::vector<mrcpp::CompFunction<3> *>> &J_ij, 
    Spinor_CompFunction &J_ij_T,
    Spinor_CompFunction &J_ij_B,
    Spinor_CompFunction &K_ij_T, 
    Spinor_CompFunction &K_ij_B,
    mrcpp::PoissonOperator &P,
    double &Orbital_energy,
    double &norm_diff_max) {

        CompFunction<3> Kappa_tree(MRA);
        CompFunction<3> Kappa_inverted_tree(MRA);
        std::vector<mrcpp::CompFunction<3> *> Nabla_Kappa_tree(3,new mrcpp::CompFunction<3>(MRA));

        CompFunction<3> Mean_field_potential(MRA);
        project(Mean_field_potential, zero, building_precision); // Initialize the potential to zero

        Spinor_gradients_pointer_list Spinor_gradients(4, std::vector<std::vector<mrcpp::CompFunction<3>*>>(2, std::vector<mrcpp::CompFunction<3>*>(3, nullptr)));
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 2; ++j) {
                for (int k = 0; k < 3; ++k) {
                    Spinor_gradients[i][j][k] = new mrcpp::CompFunction<3>(MRA);
                }
            }
        }

        std::vector<double> Energy_eigenvalues(4, 0.0);
        std::vector<std::vector<std::complex<double>>> F_matrix(4, std::vector<std::complex<double>>(4, std::complex<double>(0.0, 0.0)));
        std::vector<double> norm_difference_list(4, 0.0);

        print_memory_usage();
        std::cout << "-> Orthonormlizing Spin Orbitals" << '\n';
        orthonormalize_orbitals(Spin_orbitals, MRA);
        print_memory_usage();
        // Initialize the variables for the SCF cycle
        std::cout << "-> Calculating the Exchange contributions" << '\n';
        print_memory_usage();
        K_Psi(K_ij_T,K_ij_B, Spin_orbitals, MRA, P);
        std::cout << "-> Calculating the Coulomb contributions" << '\n';
        print_memory_usage();
        J_Psi(J_ij, Spin_orbitals, Mean_field_potential,  MRA, P);
        std::cout << "-> Calculating the gradients of the Spinors" << '\n';
        Update_Nabla_Psi_2c(Spinor_gradients, Spin_orbitals, MRA);
        print_memory_usage();
        std::cout << "-> Calculating the Kappa tree and its gradient" << '\n';
        Update_Kappa(Kappa_tree, Kappa_inverted_tree, Nabla_Kappa_tree,MRA, core_el_tree, Mean_field_potential);
        print_memory_usage();
    

        std::cout << "-> Updating the Fock matrix" << '\n';
        Update_Fock(Energy_eigenvalues, F_matrix, MRA, Spinor_gradients, Kappa_tree, Nabla_Kappa_tree, Spin_orbitals, core_el_tree, J_ij_T, J_ij_B, K_ij_T, K_ij_B);
        print_memory_usage();
        // Display the energy
        std::cout << '\t' <<  "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << '\n';
        Orbital_energy = 0.0;
        Orbital_energy = std::accumulate(Energy_eigenvalues.begin(), Energy_eigenvalues.end(), 0.0);
        std::cout << '\t' <<  "Total energy of the system = " << Orbital_energy << '\n';
        std::cout << '\t' <<  "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << '\n';
        std::cout << '\n' << '\n';
        std::cout << "Energies of the orbitals:" << '\n';
        for (const auto &energy : Energy_eigenvalues) {
            std::cout << energy << '\n';
        }

        
        std::cout << "-> Applying the Helmholtz operator to the Spinors" << '\n';
        apply_Helmholtz_ZORA_all_electrons(MRA, F_matrix, Spin_orbitals, Spinor_gradients, core_el_tree, J_ij_T, J_ij_B, K_ij_T, K_ij_B, Nabla_Kappa_tree, Kappa_tree, Kappa_inverted_tree, norm_difference_list);
        print_memory_usage();
        // Retrieve the maximum element of the norm_difference_list
        norm_diff_max = *std::max_element(norm_difference_list.begin(), norm_difference_list.end());
        

        //Delete the Spinor_gradients
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 2; ++j) {
                for (int k = 0; k < 3; ++k) {
                    Spinor_gradients[i][j][k]->free();
                    delete Spinor_gradients[i][j][k];
                }
            }
        }
        // Delete the Nabla_Kappa_tree
        for (auto &nabla_kappa : Nabla_Kappa_tree) {
            nabla_kappa->free();
            delete nabla_kappa;
        }




}




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
    Z = 4;
    n_electrons = 4; // --> As we are dealing with the H type atom
    

    std::cout << "Parameters read from file: " << '\n' << '\n';
    std::cout << " Order =" << '\t'<< '\t' << order << '\n';
    std::cout << " MaxLevel =" << '\t'<< '\t' << MaxLevel << '\n';
    std::cout << " Building precision =" << '\t' << building_precision << '\n';
    std::cout << " Epsilon =" << '\t'<< '\t' << epsilon << '\n';
    std::cout << " Z =" << '\t'<< '\t' << '\t' << Z << '\n';
    std::cout << " Number of electrons =" << '\t' << n_electrons << '\n';
    
    
    
    
    
    
    
    
    
    
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
    mrcpp::CompFunction<3> Trial_function_tree_top(mra); // -> This is the tree that will hold the output of the convolution
    mrcpp::CompFunction<3> Trial_function_tree_bottom(mra); // -> This is the tree that will hold the output of the convolution
    mrcpp::CompFunction<3> core_el_tree(mra); // -> This is the tree that will hold the [nucleus-electron] potential
    mrcpp::CompFunction<3> J_function_tree(mra);  // -> This is the tree that will hold <b a | b a> = <\psi | J - K | \psi>



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
    mrcpp::Coord<3> pos = {0,0,0};
    std::array<int,3> pow = {0,0,0};
    
    // Define the trial function as a Slater-type orbital
    std::function<double(const Coord<3> &x)> Be_1s = [] (const mrcpp::Coord<3> &r) -> double {
        auto R = std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        double alpha = 4.0;
        return std::exp(-alpha*R);
        //return std::exp(-alpha*R*R);
    };

    std::function<double(const Coord<3> &x)> Be_2s = [] (const mrcpp::Coord<3> &r) -> double {
        auto R = std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        double alpha = 1.7;
        double pref = 1-alpha*R;
        return  pref * std::exp(-0.5 *alpha*R);
        //return   std::exp(-0.5 *alpha*R*R);
    };

    


    

    // 4) Let's DEFINE now the POTENTIAL for CORE ELECTRON function V_ee(r) = -1/r
    std::function<double(const Coord<3> &x)> V_ce = [] (const mrcpp::Coord<3> &r) -> double {
        auto R = std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        return -Z/R;
    };

    // Here is the smeared potential
    // Define the smeared potential
    std::function<double(const Coord<3> &x)> smeared_potential = [] (const mrcpp::Coord<3> &r) -> double {
        auto R = std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        return -SmearedPotential::coulomb_HFYGB(R, Z, building_precision);
    };



    // We now project the potential on the tree
    //mrcpp::project(core_el_tree, V_ce, building_precision);
    mrcpp::project(core_el_tree, V_ce, building_precision);


    // 5) PROJECT the function ON the TREE
    std::vector<mrcpp::CompFunction<3>> Psi_1(2);
    std::vector<mrcpp::CompFunction<3>> Psi_2(2);
    std::vector<mrcpp::CompFunction<3>> Psi_3(2);
    std::vector<mrcpp::CompFunction<3>> Psi_4(2);



    // INITIALIZING TRIAL FUNCTIONS
    // 1s spin up
    mrcpp::project( Psi_1[0], Be_1s, building_precision);
    mrcpp::project( Psi_1[1], zero, building_precision);
    //std::cout << "Psi_2c_1[0] norm = " << Psi_1[0].getSquareNorm() << '\n';
    //std::cout << "Psi_2c_1[1] norm = " << Psi_1[1].getSquareNorm() << '\n';

    // 1s spin down
    project(Psi_2[0], zero, building_precision);
    project(Psi_2[1], Be_1s, building_precision);
    //std::cout << "Psi_2c_2[0] norm = " << Psi_2[0].getSquareNorm() << '\n';
    //std::cout << "Psi_2c_2[1] norm = " << Psi_2[1].getSquareNorm() << '\n';
    // 2s spin up
    project(Psi_3[0], Be_2s, building_precision);
    project(Psi_3[1], zero, building_precision);
    //std::cout << "Psi_2c_3[0] norm = " << Psi_3[0].getSquareNorm() << '\n';
    //std::cout << "Psi_2c_3[1] norm = " << Psi_3[1].getSquareNorm() << '\n';
    // 2s spin down
    project(Psi_4[0], zero, building_precision);
    project(Psi_4[1], Be_2s, building_precision);
    //std::cout << "Psi_2c_4[0] norm = " << Psi_4[0].getSquareNorm() << '\n';
    //std::cout << "Psi_2c_4[1] norm = " << Psi_4[1].getSquareNorm() << '\n';

    // This object is fondamental

    Orbital_Pointer_List Spin_orbitals = {
        &Psi_1,
        &Psi_2,
        &Psi_3,
        &Psi_4
    };




    


    
/*     orthonormalize_orbitals(Spin_orbitals,mra);

    

    // Now we normalize all spin orbitals
    for (auto &orbital : Spin_orbitals) {
        Renormalize_Spinor(*orbital);
        std::cout << "Spin orbital top norm =" << (*orbital)[0].getSquareNorm() << '\n';
        std::cout << "Spin orbital bottom norm =" << (*orbital)[1].getSquareNorm() << '\n';
        std::cout << '\n';
    }
   
 */  
    

    




    
    //std::vector<mrcpp::FunctionTree<3,double>> J_ij;
    std::vector<std::vector<mrcpp::CompFunction<3>*>> J_ij(4, std::vector<mrcpp::CompFunction<3>*>(2, nullptr));
    std::vector<std::vector<mrcpp::CompFunction<3>*>> K_ij(4, std::vector<mrcpp::CompFunction<3>*>(2, nullptr));
    std::vector<mrcpp::CompFunction<3>> J_ij_T(4);
    std::vector<mrcpp::CompFunction<3>> J_ij_B(4);
    std::vector<mrcpp::CompFunction<3>> K_ij_T(4);
    std::vector<mrcpp::CompFunction<3>> K_ij_B(4);
    
    for (int k =0; k<4; k++){
        project(J_ij_T[k], zero, building_precision);
        project(J_ij_B[k], zero, building_precision);
        project(K_ij_T[k], zero, building_precision);
        project(K_ij_B[k], zero, building_precision);
        J_ij[k][0] = &J_ij_T[k];
        J_ij[k][1] = &J_ij_B[k];
        K_ij[k][0] = &K_ij_T[k];
        K_ij[k][1] = &K_ij_B[k];   


        //std::cout << "J_ij[" << k << "][0] address = " << J_ij[k][0]->CompD[0] << '\n';
        //std::cout << "J_ij[" << k << "][1] address = " << J_ij[k][1]->CompD[0] << '\n';
        //std::cout << "K_ij[" << k << "][0] address = " << K_ij[k][0]->CompD[0] << '\n';
        //std::cout << "K_ij[" << k << "][1] address = " << K_ij[k][1]->CompD[0] << '\n';

    }
    





    // Create a 4x4 matrix of ComplexDouble and initialize it as zero
    std::vector<std::vector<std::complex<double>>> F_matrix(4, std::vector<std::complex<double>>(4, std::complex<double>(0.0, 0.0)));
    std::vector<double> Energy_eigenvalues(4, 0.0);
    

    Spinor_gradients_pointer_list Spinor_gradients(4, std::vector<std::vector<mrcpp::CompFunction<3>*>>(2, std::vector<mrcpp::CompFunction<3>*>(3, nullptr)));

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 2; ++j) {
            for (int k = 0; k < 3; ++k) {
                Spinor_gradients[i][j][k] = new mrcpp::CompFunction<3>(mra);
            }
        }
    }


    

    mrcpp::CompFunction<3> Mean_field_potential(mra); // -> This is the tree that will hold the potential function
    project(Mean_field_potential, zero, building_precision); // Initialize the potential to zero
    mrcpp::CompFunction<3> Kappa_tree(mra);  // -> This is the tree that will hold the K function
    mrcpp::CompFunction<3> Kappa_inverted_tree(mra);  // -> This is the tree that will hold the K^-1 function
    std::vector<mrcpp::CompFunction<3> *> Nabla_Kappa_tree(3, new mrcpp::CompFunction<3>(mra)); // -> This is the tree that will hold the gradient of K
    
    // Gradient of K
    mrcpp::ABGVOperator<3> D(mra, 0.0, 0.0); // deine the ABGV operator
   

    PoissonOperator P(mra, building_precision);



    

    //K_Psi(K_ij_T,K_ij_B, Spin_orbitals, mra, P);
    //std::cout << "K_ij_T and K_ij_B computed." << '\n';
    //std::cout << '\n' << '\n';
    //J_Psi(J_ij, Spin_orbitals, Mean_field_potential,  mra, P);
    //std::cout << "J_ij computed." << '\n';
    //std::cout << '\n' << '\n';
    //Update_Nabla_Psi_2c(Spinor_gradients, Spin_orbitals, mra);
    //std::cout << '\n' << '\n';



    //Update_Kappa(Kappa_tree, Kappa_inverted_tree, Nabla_Kappa_tree,mra, core_el_tree, Mean_field_potential);
    //std::cout << '\n' << '\n';

//    std::cout << "kin 1s" << - 0.5 * compute_Term1_T_ZORA(mra, Spinor_gradients[0], Kappa_tree, *Spin_orbitals[0]) << '\n';
//    std::cout << "Spinor gradients computed." << '\n';



    //std::cout << "Kappa tree and Nabla Kappa tree computed." << '\n';
    //Update_Fock(Energy_eigenvalues, F_matrix, mra, Spinor_gradients, Kappa_tree, Nabla_Kappa_tree, Spin_orbitals, core_el_tree, J_ij_T, J_ij_B, K_ij_T, K_ij_B);
    //std::cout << '\n' << '\n';
    //std::cout << "Fock matrix computed." << '\n';
   

    //std::cout<< '\n';
    //std::cout<< '\n';
    //std::cout << "Fock matrix:" << '\n';
    //for (const auto &row : F_matrix) {
    //    for (const auto &elem : row) {
    //        std::cout << elem << " ";
    //    }
    //    std::cout << '\n';
    //}
    //std::cout<< '\n';
    //std::cout<< '\n';

    //std::cout << "Energy eigenvalues:" << '\n'; 
    //for (const auto &energy : Energy_eigenvalues) {
    //    std::cout << energy << '\n';
    //}
    

    //std::cout << '\n' << '\n';
    //std::cout << "************************************************************"  << '\n';
    //std::cout << "Tot energy of the system = " << std::accumulate(Energy_eigenvalues.begin(), Energy_eigenvalues.end(), 0.0) << '\n';
    //std::cout << "************************************************************" << '\n' << '\n';

    std::vector<double> norm_difference_list(4,0.0);

    //apply_Helmholtz_ZORA_all_electrons(mra, F_matrix, Spin_orbitals, Spinor_gradients, core_el_tree, J_ij_T, J_ij_B, K_ij_T, K_ij_B, Nabla_Kappa_tree, Kappa_tree, Kappa_inverted_tree, norm_difference_list);
    //std::cout << '\n' << '\n';


    //std::cout << "New psi norms:" << '\n';
    //for (int i = 0; i < 4; ++i) {
    //    std::cout << "Psi_2c_" << i << "[0] norm = " << (*Spin_orbitals[i])[0].getSquareNorm() << '\t';
    //    std::cout << "Psi_2c_" << i << "[1] norm = " << (*Spin_orbitals[i])[1].getSquareNorm() << '\n';
    //}


    //orthonormalize_orbitals(Spin_orbitals, mra);

    //std::cout << '\n' << '\n';
    //std::cout << "New round." << '\n';

    //orthonormalize_orbitals(Spin_orbitals, mra);

    
    //std::cout << "Orthonormalized orbitals." << '\n';
    //std::cout << "New psi norms after orthonormalization:" << '\n';
    //for (int i = 0; i < 4; ++i) {
    //    std::cout << "Psi_2c_" << i << "[0] norm = " << (*Spin_orbitals[i])[0].getSquareNorm() << '\t';
    //    std::cout << "Psi_2c_" << i << "[1] norm = " << (*Spin_orbitals[i])[1].getSquareNorm() << '\n';
    //}




    //K_Psi(K_ij_T,K_ij_B, Spin_orbitals, mra, P);
    //std::cout << '\n' << '\n';
    //J_Psi(J_ij, Spin_orbitals, Mean_field_potential,  mra, P);
    //std::cout << '\n' << '\n';
    //Update_Nabla_Psi_2c(Spinor_gradients, Spin_orbitals, mra);
    //std::cout << '\n' << '\n';
    //Update_Kappa(Kappa_tree, Kappa_inverted_tree, Nabla_Kappa_tree,mra, core_el_tree, Mean_field_potential);
    //std::cout << '\n' << '\n';
    //Update_Fock(Energy_eigenvalues, F_matrix, mra, Spinor_gradients, Kappa_tree, Nabla_Kappa_tree, Spin_orbitals, core_el_tree, J_ij_T, J_ij_B, K_ij_T, K_ij_B);
    //std::cout << '\n' << '\n';
   //std::cout<< '\n';
    //std::cout<< '\n';
    //std::cout << "Energy eigenvalues:" << '\n'; 
    //for (const auto &energy : Energy_eigenvalues) {
    //    std::cout << energy << '\n';
    //}
    //std::cout << '\n' << '\n';
    //std::cout << "************************************************************"  << '\n';
    //std::cout << "Tot energy of the system = " << std::accumulate(Energy_eigenvalues.begin(), Energy_eigenvalues.end(), 0.0) << '\n';
    //std::cout << "************************************************************" << '\n' << '\n';













    double Psi_norm = 0;

    


    

   
      /*
     *-----------------------------------
     *-----------------------------------
     *         BEGIN SCF CYCLE
     *-----------------------------------
     *-----------------------------------
    */

    // Parameters for the SCF
    double norm_diff_max = 1.0;   //    -> Norm of the difference in the 2 consecutive iterations (initialized as 1 to begin the while loop)
    double E_tot;               //    -> Energy of the system
    double E_Psi;               //    -> Energy of the orbital
    
    // We now define all the trees that will be used to compute the enrgy and the SCF cycle
     // Using this to calculate the gradient of K, in the one el case is a constant.
    //Nabla_K_tree = mrcpp::gradient(D, K_tree);
    std::vector<std::vector<mrcpp::CompFunction<3>*>> Nabla_Psi_2c(2, std::vector<mrcpp::CompFunction<3>*>(3, new mrcpp::CompFunction<3>(mra))); // -> This is the tree that will hold the gradient of Psi_2c

    
    // ==================================================================================================================================
    
    
    // Create a vector to hold both Nabla_Psi_t and Nabla_Psi_b
    
    
    std::cout << "************************************************************" << '\n';
    
    
    std::cout << '\n' << '\n';
    // A few utilities variables for the SCF cycle
    Orbital_Pointer_List Spin_orbitals_next = {
        nullptr,
        nullptr,
        nullptr,
        nullptr};   // -> This will hold the next iteration of the spinor
    mrcpp::CompFunction<3> Psi_2c_diff_top(mra);                // -> This will hold the difference between the 2 spinors for the TOP component
    mrcpp::CompFunction<3> Psi_2c_diff_bottom(mra);             // -> This will hold the difference between the 2 spinors for the BOTTOM component
    
    


    int max_cycle = 9;
    
    
    
    //===============================================================================

 

    std::vector<double> E_values;
    std::vector<double> norm_diff_values;
    std::vector<int> cycle_values;

    

    while (norm_diff_max > epsilon) {
        num_cycle++;    
        std::cout << "************************************************************" << '\n';
        std::cout << "> CYCLE " << num_cycle << " STARDED:" << '\n' << '\n';
        
        SCF_Cycle_ZORA( mra,
                        Spin_orbitals,
                        core_el_tree,
                        J_ij,
                        J_ij_T,
                        J_ij_B,
                        K_ij_T,
                        K_ij_B,
                        P,
                        E_tot,
                        norm_diff_max);

                    


        /*
        print_memory_usage();
        std::cout << "-> Orthonormlizing Spin Orbitals" << '\n';
        orthonormalize_orbitals(Spin_orbitals, mra);
        print_memory_usage();
        // Initialize the variables for the SCF cycle
        std::cout << "-> Calculating the Exchange contributions" << '\n';
        print_memory_usage();
        K_Psi(K_ij_T,K_ij_B, Spin_orbitals, mra, P);
        std::cout << "-> Calculating the Coulomb contributions" << '\n';
        print_memory_usage();
        J_Psi(J_ij, Spin_orbitals, Mean_field_potential,  mra, P);
        std::cout << "-> Calculating the gradients of the Spinors" << '\n';
        Update_Nabla_Psi_2c(Spinor_gradients, Spin_orbitals, mra);
        print_memory_usage();
        std::cout << "-> Calculating the Kappa tree and its gradient" << '\n';
        Update_Kappa(Kappa_tree, Kappa_inverted_tree, Nabla_Kappa_tree,mra, core_el_tree, Mean_field_potential);
        print_memory_usage();
    

        // Compute the energy
        if (verbose) {std::cout << "$ STARTING TO COMPUTE THE ENERGY..." << '\n'<< '\n';}

        std::cout << "-> Updating the Fock matrix" << '\n';
        Update_Fock(Energy_eigenvalues, F_matrix, mra, Spinor_gradients, Kappa_tree, Nabla_Kappa_tree, Spin_orbitals, core_el_tree, J_ij_T, J_ij_B, K_ij_T, K_ij_B);
        print_memory_usage();
        // Display the energy
        std::cout << '\t' <<  "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << '\n';
        E_tot = 0.0;
        E_tot = std::accumulate(Energy_eigenvalues.begin(), Energy_eigenvalues.end(), 0.0);
        std::cout << '\t' <<  "Total energy of the system = " << E_tot << '\n';
        std::cout << '\t' <<  "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << '\n';
        std::cout << '\n' << '\n';
        std::cout << "Energies of the orbitals:" << '\n';
        for (const auto &energy : Energy_eigenvalues) {
            std::cout << energy << '\n';
        }





    
        // Apply the Helmholtz operator
        if (verbose) {std::cout << "$ STARTING TO COMPUTE NEXT ITERATION'S WAVEFUNCTION..." << '\n'<< '\n';}
        std::cout << "-> Applying the Helmholtz operator to the Spinors" << '\n';
        apply_Helmholtz_ZORA_all_electrons(mra, F_matrix, Spin_orbitals, Spinor_gradients, core_el_tree, J_ij_T, J_ij_B, K_ij_T, K_ij_B, Nabla_Kappa_tree, Kappa_tree, Kappa_inverted_tree, norm_difference_list);
        print_memory_usage();
        // Retrieve the maximum element of the norm_difference_list
        norm_diff_max = *std::max_element(norm_difference_list.begin(), norm_difference_list.end());
        */
        
        std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << '\n';
        std::cout << "@ Norm of the difference = " << norm_diff_max << " @" << '\n';
        std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << '\n';
        std::cout << '\n';
        
        
        
       
        std::cout << "SCF CYCLE COMPLETED!" << '\n' << '\n';


        // Store all the values for the SCF cycle
        E_values.push_back(E_tot);
        norm_diff_values.push_back(norm_diff_max);
        cycle_values.push_back(num_cycle);

        std::cout << "Last cleanup, MOVING ON TO THE NEXT CYCLE!" << '\n' << '\n';



        if (num_cycle >= max_cycle){
            std::cout << "The SCF cycle did not converge" << '\n';
            print::footer(0, timer, 2);
            exit(0);
        }
        
}



    std::cout << "************************************************************" << '\n';


    std::cout << '\n' << '\n';
    std::cout << "-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-" << '\n';
    std::cout << "Your calculation is finished, HURRAY. Here are some cute cats:"<<'\n' << '\n';
    std::cout << "=^.^=" << '\t' <<"=^.^=" << '\t' <<"=^.^=" << '\t' <<"=^.^=" << '\t' <<"=^.^=" << '\t' <<"=^.^=" << '\t' <<"=^.^=" << '\t' << '\n' << '\n';

    std::cout << "Here is a summary of the SCF cycle:" << '\n' << '\n';
    std::cout << "Cycle" << '\t' << "Energy" << '\t' << '\t' << '\t' << "Norm. Diff." << '\n';
    for (int i = 0; i < E_values.size(); i++){
        std::cout << cycle_values[i] << '\t' << E_values[i] << '\t' << norm_diff_values[i] << '\n';
    }
    std::cout << '\n' << '\n' << "Non relativistic energuy for He atom" << '\n';
    std::cout << "*************************************************************" << '\n';
    //(mra, *He_NR_Solution.CompD[0], *core_el_tree.CompD[0], mu, E_tot);
    std::cout << "*************************************************************" << '\n';

    print::footer(0, timer, 2);
    return 0;
}