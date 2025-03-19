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














// ==================================================================================================================================


double Energy_Non_Relativistic(double precision, MultiResolutionAnalysis<3> &MRA, FunctionTree<3> &f, FunctionTree<3> &V_ce, FunctionTree<3> &J, int &n_electrons){
    FunctionTree<3> f_V(MRA);
    FunctionTreeVector<3> Nabla_f; // THE GRADIENT WILL BE A VECTOR!!!!!!!
    double norm_sqrd = dot(f,f);
 
    mrcpp::ABGVOperator<3> D(MRA, 0.0, 0.0);
    
    // ONE ELECTRON PART
    Nabla_f = mrcpp::gradient(D,f);

    FunctionTree<3> Nabla_fx(MRA); 
    FunctionTree<3> Nabla_fy(MRA); 
    FunctionTree<3> Nabla_fz(MRA); 
    get_func(Nabla_f,0).deep_copy(&Nabla_fx);
    get_func(Nabla_f,1).deep_copy(&Nabla_fy);
    get_func(Nabla_f,2).deep_copy(&Nabla_fz);

    // Compute the kinetic energy
    double nabla_f_norm_sqrd = dot(Nabla_fx,Nabla_fx) + dot(Nabla_fy,Nabla_fy) + dot(Nabla_fz,Nabla_fz);
    double kinetic_energy = 0.5000000000 * nabla_f_norm_sqrd/ norm_sqrd;
    
    // Compute the product of the gradient and the function
    mrcpp::multiply(precision,f_V, 1.0, f, V_ce);
    double potential_energy_one_electron = dot(f,f_V) / norm_sqrd;

    double One_electron_energy = potential_energy_one_electron + kinetic_energy;

    // TWO ELECTRON PART - Only applyable if we have more than 2 electrons
    double Two_electron_energy = 0;
    if (n_electrons == 2){
        FunctionTree<3> Jf_tree(MRA);
        mrcpp::multiply(precision,Jf_tree, 1.0, J, f);
        // Two_electron_energy = \sum_{ij} < \Phi_i | J | \Phi_j > - < \Phi_i | K | \Phi_j > = < f | Jf > - 0 
        Two_electron_energy = n_electrons * dot(f,Jf_tree) / norm_sqrd; // should be for each of the 2 electrons but they occupy the same spatial orbital, being orthonormal for the spin
        // multiplied by 2 because we have 2 electrons, so one for each spin-orbital
    }


    double potential_energy_total = potential_energy_one_electron * n_electrons + 0.5 * Two_electron_energy;

    std::cout << '\n';
    std::cout << " T    = " << n_electrons * kinetic_energy << '\n';
    std::cout << " V_ce = " << n_electrons * potential_energy_one_electron  << '\n';
    std::cout << " V_ee = " << 0.5 * Two_electron_energy <<'\n';
    std::cout << " V_t  = " << potential_energy_total << '\n' << '\n';

    double result = One_electron_energy + 0.5 * Two_electron_energy;

    std::cout << "  ########################## " << '\n';
    std::cout << "  | E = " << result << " |" << '\n';
    std::cout << "  ########################## " << '\n' << '\n';
    return result;
}

// ==================================================================================================================================


double Energy_ZORA(double precision, MultiResolutionAnalysis<3> &MRA, FunctionTree<3> &f, FunctionTree<3> &V){
    // Declare the FunctionTrees
    FunctionTree<3> f_V(MRA);
    FunctionTree<3> K_tree(MRA);
    FunctionTree<3> Kf_tree(MRA);

    // Define the K function
    std::function<double(const Coord<3> &x)> K_r = [&V] (const mrcpp::Coord<3> &r) -> double {
        double result = 1.0 - V.evalf(r) / (2.0 * m * c * c);
        return 1.0 / result;
    };

    // And project it on the correspondig tree
    mrcpp::project<3>(precision, K_tree, K_r);



    // Define the type of derivative operator: ABGV
    mrcpp::ABGVOperator<3> D(MRA, 0.0, 0.0);


    // *********************************************
    // Now we can start the subroutine:
    // Compute the norm of the function (not strictly needed)
    double norm_sqrd = dot(f,f);
 

    // Compute the gradient of f
    // Now we declare all the function trees needed
    FunctionTreeVector<3> Nabla_f; 
    FunctionTreeVector<3> Nabla_K;
    FunctionTreeVector<3> Nabla_Kf;
    Nabla_f = mrcpp::gradient(D,f);
    Nabla_K = mrcpp::gradient(D,K_tree);
    
    // Compute K(r) * f
    mrcpp::multiply(precision,Kf_tree, 1.0, K_tree, f);

    // And finde the gradient of it
    Nabla_Kf = mrcpp::gradient(D,Kf_tree);

    // Assign the components to the corresponding functions
    FunctionTree<3> Nabla_fx(MRA); 
    FunctionTree<3> Nabla_fy(MRA); 
    FunctionTree<3> Nabla_fz(MRA); 
    FunctionTree<3> Nabla_Kfx(MRA);
    FunctionTree<3> Nabla_Kfy(MRA);   
    FunctionTree<3> Nabla_Kfz(MRA);
    get_func(Nabla_f,0).deep_copy(&Nabla_fx);
    get_func(Nabla_f,1).deep_copy(&Nabla_fy);
    get_func(Nabla_f,2).deep_copy(&Nabla_fz);
    get_func(Nabla_Kf,0).deep_copy(&Nabla_Kfx);
    get_func(Nabla_Kf,1).deep_copy(&Nabla_Kfy);
    get_func(Nabla_Kf,2).deep_copy(&Nabla_Kfz);


    // Same for the K gradient
    FunctionTree<3> Nabla_Kx(MRA);
    FunctionTree<3> Nabla_Ky(MRA);
    FunctionTree<3> Nabla_Kz(MRA);
    get_func(Nabla_K,0).deep_copy(&Nabla_Kx);
    get_func(Nabla_K,1).deep_copy(&Nabla_Ky);
    get_func(Nabla_K,2).deep_copy(&Nabla_Kz);

    // p \cdot K(r) p = - < f \Nabla(K) | \Nabla(f) > + < \Nabla(f * K) | \Nabla(f) >

    // Compute the first term of the ZORA kinetic energy, but first we need the [f * \Nabla(K)]
    FunctionTree<3> f_Nabla_Kx(MRA);
    FunctionTree<3> f_Nabla_Ky(MRA);
    FunctionTree<3> f_Nabla_Kz(MRA);
    mrcpp::multiply(precision,f_Nabla_Kx, 1.0, f, Nabla_Kx);
    mrcpp::multiply(precision,f_Nabla_Ky, 1.0, f, Nabla_Ky);
    mrcpp::multiply(precision,f_Nabla_Kz, 1.0, f, Nabla_Kz);

    double f_Nabla_K__dot__Nabla_f = dot(f_Nabla_Kx,Nabla_fx) + dot(f_Nabla_Ky,Nabla_fy) + dot(f_Nabla_Kz,Nabla_fz);

    // Compute the second term of the ZORA kinetic energy
    double Nabla_Kf__dot__Nabla_f = dot(Nabla_Kfx,Nabla_fx) + dot(Nabla_Kfy,Nabla_fy) + dot(Nabla_Kfz,Nabla_fz);

    
    // Compute the kinetic energy
    double kinetic_energy = 1 / (2 * m);
    kinetic_energy = kinetic_energy * ( -f_Nabla_K__dot__Nabla_f + Nabla_Kf__dot__Nabla_f);
    


    // For the potential part:
    mrcpp::multiply(precision,f_V, 1.0, f, V);
    double potential_energy = dot(f,f_V) / norm_sqrd;

    double result = potential_energy + kinetic_energy;

    std::cout << "T_ZORA = " << kinetic_energy << '\n';
    std::cout << "V = " << potential_energy << '\n' << '\n';

    std::cout << "  ########################### " << '\n';
    std::cout << "  | E = " << result << " |" << '\n';
    std::cout << "  ########################### " << '\n' << '\n';
    return result;
}

// ==================================================================================================================================

double energy_update(MultiResolutionAnalysis<3> &MRA, FunctionTree<3> &Tilde_phi, FunctionTree<3> &Phi, FunctionTree<3> &V){
    FunctionTree<3> DeltaPhi(MRA);  // Allocate space for the difference
    FunctionTree<3> V_DeltaPhi(MRA); // Allocate space for the product of the potential and the difference

    // norm = <Tilde_phi|Tilde_phi>
    double norm = dot(Tilde_phi,Tilde_phi);
 
    // |DeltaPhi> = |Tilde_phi> - |Phi>
    mrcpp::add(1e-5,DeltaPhi,-1.0,Phi, 1.0,Tilde_phi);
    
    // |V_DeltaPhi> = V * |DeltaPhi>
    mrcpp::multiply(0.000001,V_DeltaPhi, 1.0, V, DeltaPhi);
    
    // <Tilde_phi|V * DeltaPhi>
    double result = dot(V_DeltaPhi,Tilde_phi);
    // <Tilde_phi|V * DeltaPhi> / <Tilde_phi|Tilde_phi>
    return (result / norm) ;
}


// ==================================================================================================================================



void apply_Helmholtz(double mu,double building_precision, MultiResolutionAnalysis<3> &MRA, FunctionTree<3> &out, FunctionTree<3> &f){
    HelmholtzOperator Helm(MRA, mu, building_precision);
    apply(building_precision, out, Helm, f, -1, false);
}



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






void Build_J_operator(mrcpp::MultiResolutionAnalysis<3> &MRA,mrcpp::FunctionTree<3> &J, mrcpp::FunctionTree<3> &Phi, int n_electrons){
    if (n_electrons == 1){
        return;
    }

    mrcpp::FunctionTree<3> rho(MRA); // -> This is the tree that will hold the electron density
    mrcpp::FunctionTree<3> J_tmp(MRA);
    // rho = Phi^2
    mrcpp::multiply(building_precision*0.01, rho, 1.0, Phi, Phi);


    // I build the Poisson operator as a L.C. of Gaussians (so it's faster to integrate)
    mrcpp::PoissonOperator Poisson(MRA, building_precision);
    J.clear();
    // J = Poisson(rho)
    mrcpp::apply(building_precision, J_tmp, Poisson, rho, -1, false);

    // J = J_tmp + J_tmp
    mrcpp::add(building_precision, J, 1.0, J_tmp, -0.5, J_tmp);
}
//==================================================================================================================================


void compute_2e_energy(MultiResolutionAnalysis<3> &MRA,int n_electrons, FunctionTreeVector<3,double> Phi_vec, FunctionTree<3> &V, FunctionTreeVector<3,double> &J_vec, FunctionTreeVector<3,double> &K_Phi_vec){
    FunctionTree<3> J_tmp(MRA);
    FunctionTree<3> J_accumulation(MRA);
    FunctionTree<3> J_tmp_sum(MRA);
    FunctionTree<3> K_tmp(MRA);
    FunctionTree<3> K_tmp_Phi_b(MRA);
    FunctionTree<3> K_accumulation(MRA);
    FunctionTree<3> K_tmp_sum(MRA);
    FunctionTree<3> Phi_tmp(MRA);
    FunctionTree<3> PhiA_star_times_PhiA(MRA);
    FunctionTree<3> PhiB_star_times_PhiB(MRA);
    FunctionTree<3> PhiA_star_times_PhiB(MRA);
    
    mrcpp::PoissonOperator Poisson(MRA, building_precision);

    std::vector<double> J_n;
    std::vector<double> K_n;

    J_n.clear();
    K_n.clear();

    // Compute the J operator for each of the electrons
    for (int a = 0; a < n_electrons/2; a++){
        PhiA_star_times_PhiA.clear();
        // \rho_a = \Phi[a]^* \Phi[a]
        
        mrcpp::multiply(building_precision,PhiA_star_times_PhiA, 1.0, get_func(Phi_vec, a), get_func(Phi_vec, a));
        
        J_accumulation.clear();
        for (int b = 0; b < n_electrons; b++){
            J_tmp_sum.clear();
            PhiB_star_times_PhiB.clear();
            // \rho_b = \Phi[b]^* \Phi[b]
            mrcpp::multiply(building_precision,PhiB_star_times_PhiB, 1.0, get_func(Phi_vec, b), get_func(Phi_vec, b));
            
            J_tmp.clear();
            // J_tmp(r) = \int \rho_b(r') / |r-r'| dr'
            mrcpp::apply(building_precision, J_tmp, Poisson, PhiB_star_times_PhiB, -1, false);    
            
            // J_tmp_sum = J_accumulation + J_tmp: J_accumulation = \sum_{i=0}^{b} J_tmp[i]
            mrcpp::add(building_precision, J_tmp_sum, 1.0, J_accumulation, 1.0, J_tmp);
            
            J_tmp_sum.deep_copy(&J_accumulation);

        }
        // save the operator for this orbital
        J_accumulation.deep_copy(&get_func(J_vec, a)); // !!!!!!! -----> This is the J operator for the a-th as just \hat{J}
        // Now we exit the first loop and we have accumulated the J operator for the a-th electron, by summing the densitie of all the electrons.

        // Expectation value for the a-th orbital:
        // J_n = < \Phi[a] | J | \Phi[a] >
        J_n.push_back(dot(PhiA_star_times_PhiA, J_accumulation));
    }


    // Compute the K operator for each of the electrons
    for (int a = 0; a < n_electrons/2; a++){
        K_accumulation.clear();
        for (int b = 0; b < n_electrons; b++){
                K_tmp_sum.clear();


                PhiA_star_times_PhiB.clear();
                // \rho_ab = \Phi[a]^* \Phi[b]
                mrcpp::multiply(building_precision,PhiA_star_times_PhiB, 1.0, get_func(Phi_vec, a), get_func(Phi_vec, b));


                K_tmp.clear(); 
                // K_tmp(r) = \int \rho_b(r') / |r-r'| dr'
                mrcpp::apply(building_precision, K_tmp, Poisson, PhiA_star_times_PhiB, -1, false);    
                
                K_tmp_Phi_b.clear();
                // K_tmp_Phi = \Phi[b] * \int \Phi[a]^* * \Phi[b] / |r-r'| dr'
                mrcpp::multiply(building_precision, K_tmp_Phi_b, 1.0, K_tmp, get_func(Phi_vec,b));


                // K_tmp_sum = K_accumulation + K_tmp_Phi:      K_accumulation = \sum_{i=0}^{b} K_tmp_Phi_b[i]
                mrcpp::add(building_precision, K_tmp_sum, 1.0, K_accumulation, 1.0, K_tmp_Phi_b);
                
                K_tmp_sum.deep_copy(&K_accumulation);
        }
        K_accumulation.deep_copy(&get_func(K_Phi_vec, a)); // !!!!!-----> This is the K operator for the a-th as \hat{K}\ket{\Phi[a]}
        // Now we exit the first loop and we have accumulated the K operator for the a-th electron

        // Expectation value for the a-th orbital:
        // K_n = < \Phi[a] | K | \Phi[a] >
        J_n.push_back(dot(get_func(Phi_vec,a), K_accumulation));
    }
    
    double Two_electron_energy;
    double coulomb_energy = std::accumulate(J_n.begin(), J_n.end(), 0.0);
    double exchange_energy = std::accumulate(K_n.begin(), K_n.end(), 0.0);
    Two_electron_energy = 2 * coulomb_energy - exchange_energy;
    std::cout << "_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-" << '\n';
    std::cout << "Computing the two electron energy contributions..." << '\n' << '\n';

    std::cout << "Coulomb energy = " << coulomb_energy << '\n';
    std::cout << "Exchange energy = " << exchange_energy << '\n';
    std::cout << "Two electron energy = " << Two_electron_energy << '\n';
    std::cout << "************************************************************" << '\n';




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
    mrcpp::FunctionTree<3> GVPhi_tree(mra); // -> This is the tree that will hold the output of the convolution
    mrcpp::FunctionTree<3> J_tree(mra); // -> This is the tree that will hold the [electron-electron] (pair) potential
    mrcpp::FunctionTree<3> core_el_tree(mra); // -> This is the tree that will hold the [nucleus-electron] potential

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
    mrcpp::GaussFunc<3> Phi_trial(beta, alpha, pos, pow);



    // 4) Let's DEFINE now the POTENTIAL for CORE ELECTRON function V_ee(r) = -1/r
    std::function<double(const Coord<3> &x)> V_ce = [] (const mrcpp::Coord<3> &r) -> double {
        auto R = std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        return -Z/R;
    };

    // 5) PROJECT the function ON the TREE
    mrcpp::project<3>(building_precision, GVPhi_tree, Phi_trial); // I project the trial-function on the tree
    GVPhi_tree.normalize();
    std::cout << "Psi_trial" << GVPhi_tree << '\n';

    // As well as the potential:
    if (n_electrons == 1){
    mrcpp::project<3>(building_precision, core_el_tree, V_ce); // I do the same with the potential
    core_el_tree.deep_copy(&Potential_tree);
    }
    else if (n_electrons == 2){
        mrcpp::project<3>(building_precision, core_el_tree, V_ce); // I do the same with the potential
        Build_J_operator(mra,J_tree, GVPhi_tree, n_electrons);
        mrcpp::add(building_precision, Potential_tree, double(n_electrons), core_el_tree, 1.0, J_tree);
    }
    else{
        std::cout << "The number of elecrons is not supported" << '\n';
        return 0;
    }



    if (debug){
    // Befofe starting the SCF cycle, let's print some debugging information:
    std::cout << "Here is some debugging information: " << '\n';
    std::cout << "************************************" << '\n';
    std::cout << "Potential_tree = " << '\n';
    std::cout << Potential_tree << '\n';
    std::cout << "************************************" << '\n';
    std::cout << "Gauss_tree = " << '\n';
    std::cout << GVPhi_tree << '\n';
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


    // Print the norm of the difference
    std::cout << "Cycle " << num_cycle << " done...  Norm of the difference = "<< norm_diff << '\n';
    if (Relativity == 0){
        E = Energy_Non_Relativistic(building_precision, mra, GVPhi_tree, core_el_tree, J_tree, n_electrons);
    }
    else if (Relativity == 1){
        E = Energy_ZORA(building_precision, mra, GVPhi_tree, Potential_tree);
    }
    else{
        std::cout << "Relativity not supported" << '\n';
        return 0;
    }
    
    //std::cout << "Energy = " << E << '\n';
    mu = std::sqrt(-2.0 * E);
    std::cout << "mu = " << mu << '\n';
    std::cout << "************************************************************" << '\n';
    num_cycle++;

    


    return 0;


    while (norm_diff > epsilon) {
        


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
        
        apply_Helmholtz(mu,building_precision, mra, GVPhi_tree, VPhi_tree);

        //std::cout << "Energy update: " << energy_update(mra, GVPhi_tree, Phi_n_copy_tree, Potential_tree) << '\n';
        

        // Normalize the function: \Tilde{\Phi^{n+1}} --> \Phi^{n+1}
        GVPhi_tree.normalize();
        Build_J_operator(mra,J_tree, GVPhi_tree, n_electrons);
        // Normalized_difference_tree = \Phi^{n+1} - \Phi^{n}
        mrcpp::add(apply_precision, Normalized_difference_tree, +1.0, GVPhi_tree, -1.0, Phi_n_copy_tree);
        // norm_diff = || \Phi^{n+1} - \Phi^{n} ||
        norm_diff = Normalized_difference_tree.getSquareNorm();
        norm_diff = std::sqrt(norm_diff);





        // Print the norm of the difference
        std::cout << "Cycle " << num_cycle << " done...  Norm of the difference = "<< norm_diff << '\n';
        if (Relativity == 0){
            E = Energy_Non_Relativistic(building_precision, mra, GVPhi_tree, core_el_tree, J_tree, n_electrons);
        }
        else if (Relativity == 1){
            E = Energy_ZORA(building_precision, mra, GVPhi_tree, Potential_tree);
        }
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

    print::footer(0, timer, 2);
    return 0;
}