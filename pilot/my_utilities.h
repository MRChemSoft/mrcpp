#pragma once

#include <cstdlib>
#include <numeric>
#include <iostream>
#include <vector>

#include "MRCPP/Printer"
#include "MRCPP/Gaussians"
#include "MRCPP/MWFunctions"
#include "MRCPP/MWOperators"


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

/* 
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

 */// ==================================================================================================================================

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




