//void compute_V_phi_with_spin_orbit_coupling(){
//}
#pragma once

#include <iostream>
#include <vector>
#include <cstdlib>
#include <numeric>


#include "../api/Printer"
#include "../api/Timer"
#include "../api/Gaussians"
#include "../api/Plotter"
#include "../api/MWFunctions"
#include "../api/MWOperators"






// FOR NOW I PLAY IT SAFE, A SIMPLE ADDITION MAY BE JUST TO ADD IN THE ARGUMENTS THE ABGV OPERATOR, NOT TO REDEFINE IT EVERY TIME



/*
*
*   CHANGE IN TERM 1 AND 2: THE FUNCTION DOT ALREADY CONJUGATES THE COMPFUNCTION BRA!!!!!!!!! UNCONJUGATE THEM AS IT IS DONE IN THE FUNCTION
*
*/


// IN principle, one could skip the step of recomputing all the gradient of the \Psi * K bu simpli adding \Nabla \Psi  K + \Psi  \Nabla K, but this would be longer, tho requiring less memory
ComplexDouble compute_Term1_T_ZORA(MultiResolutionAnalysis<3> &MRA, std::vector<std::vector<mrcpp::CompFunction<3>*>> &Nabla_Psi_2c, mrcpp::CompFunction<3> &K_tree,  std::vector<mrcpp::CompFunction<3>> Psi_2c){
    // Nabla(\Psi K) = \Nabla \Psi  K + \Psi  \Nabla K

    CompFunction<3> Psi_t_K;
    CompFunction<3> Psi_b_K;

    std::vector<mrcpp::CompFunction<3> *> Nabla_Psi_t = Nabla_Psi_2c[0];
    std::vector<mrcpp::CompFunction<3> *> Nabla_Psi_b = Nabla_Psi_2c[1];

    // Compute Psi_top * K
    mrcpp::multiply(Psi_t_K, Psi_2c[0], K_tree,building_precision, false, false, false);
    // Compute Psi_bottom * K
    mrcpp::multiply(Psi_b_K, Psi_2c[1], K_tree,building_precision, false, false, false);

    // Operatr ABGV
    mrcpp::ABGVOperator<3> D(MRA, 0.0, 0.0);

    // Gradient of Psi_top * K
    auto Nabla_Psi_t_K = mrcpp::gradient(D,Psi_t_K);
    // Gradient of Psi_bottom * K
    auto Nabla_Psi_b_K = mrcpp::gradient(D,Psi_b_K);

    // REMEMBER THE MINUS SIGN!!! (it comes from the integration by parts)
    ComplexDouble Top_contribute = -1.0*(dot(*Nabla_Psi_t_K[0],*Nabla_Psi_t[0]) + dot(*Nabla_Psi_t_K[1],*Nabla_Psi_t[1]) + dot(*Nabla_Psi_t_K[2],*Nabla_Psi_t[2]));
    ComplexDouble Bottom_contribute = -1.0*(dot(*Nabla_Psi_b_K[0],*Nabla_Psi_b[0]) + dot(*Nabla_Psi_b_K[1],*Nabla_Psi_b[1]) + dot(*Nabla_Psi_b_K[2],*Nabla_Psi_b[2]));
    return Top_contribute + Bottom_contribute;
}



ComplexDouble compute_Term2_T_ZORA(MultiResolutionAnalysis<3> &MRA, std::vector<std::vector<mrcpp::CompFunction<3>*>> &Nabla_Psi_2c, mrcpp::CompFunction<3> &K_tree, std::vector<mrcpp::CompFunction<3> *> &Nabla_K_tree,  std::vector<mrcpp::CompFunction<3>> Psi_2c){
    std::vector<mrcpp::CompFunction<3>> Psi_t_Nabla_K;; 
    std::vector<mrcpp::CompFunction<3>> Psi_b_Nabla_K;

    CompFunction<3> Psi_Nabla_k_tmp(MRA);
    CompFunction<3> Psi_element(MRA);
    CompFunction<3> Nabla_K_element(MRA);
    
    // Operatr ABGV
    mrcpp::ABGVOperator<3> D(MRA, 0.0, 0.0);


    // Gradient of K
    auto Nabla_K = Nabla_K_tree;

    // Gradient of the 2 components of Psi
    // PAY ATTENTION!!! The gradient obtained is actually a pointer, so later on we'll have to dereference it
    auto Nabla_Psi_t = Nabla_Psi_2c[0];
    auto Nabla_Psi_b = Nabla_Psi_2c[1];

    // Compute Psi_top * Nabla_K
    for (int i=0;i<3;i++){
        mrcpp::multiply(Psi_Nabla_k_tmp, Psi_2c[0], *Nabla_K[i],building_precision, false, false, false);   // <---- Here we dereference the pointer
        Psi_t_Nabla_K.push_back(Psi_Nabla_k_tmp);
    }
    
    // Compute Psi_bottom * Nabla_K
    for (int i=0;i<3;i++){
        mrcpp::multiply(Psi_Nabla_k_tmp, Psi_2c[1], *Nabla_K[i],building_precision, false, false, true);   // <---- Here we dereference the pointer
        Psi_b_Nabla_K.push_back(Psi_Nabla_k_tmp);
    }

    std::vector<ComplexDouble> Kinetic_term_1;


    ComplexDouble Top_contribution  = dot(Psi_t_Nabla_K[0],*Nabla_Psi_t[0]) + dot(Psi_t_Nabla_K[1],*Nabla_Psi_t[1]) + dot(Psi_t_Nabla_K[2],*Nabla_Psi_t[2]);
    ComplexDouble Bottom_contribution  = dot(Psi_b_Nabla_K[0],*Nabla_Psi_b[0]) + dot(Psi_b_Nabla_K[1],*Nabla_Psi_b[1]) + dot(Psi_b_Nabla_K[2],*Nabla_Psi_b[2]);


    return Top_contribution + Bottom_contribution;
}









void compute_rotor( std::vector<mrcpp::CompFunction<3>*> &Rotor , mrcpp::ABGVOperator<3> &D,  std::vector<mrcpp::CompFunction<3>*> &K_Nabla_Psi, MultiResolutionAnalysis<3> &MRA){
    Eigen::Matrix<int, 3, 2> Curl_coef;
    Curl_coef << 1, 2,
                 2, 0,
                 0, 1;
                      
    CompFunction<3> Deriv_1(MRA);
    CompFunction<3> Deriv_2(MRA);


    // For synthax:
    //void mrcpp::apply<3>(mrcpp::CompFunction<3> &out, mrcpp::DerivativeOperator<3> &oper, mrcpp::CompFunction<3> &inp, int dir, ComplexDouble (*metric)[4] = (ComplexDouble (*)[4])nullptr)
    // mrcpp::apply(deriv_Psi1, D, Psi_in[s], j); // Derivative of Psi with respect to j
    for (int i=0; i<3; i++){
        mrcpp::apply<3>(Deriv_1, D, *(K_Nabla_Psi[Curl_coef(i,1)]), Curl_coef(i,0));
        mrcpp::apply<3>(Deriv_2, D, *(K_Nabla_Psi[Curl_coef(i,0)]), Curl_coef(i,1));
        mrcpp::add<3>(*Rotor[i], 1.0, Deriv_1, -1.0, Deriv_2, building_precision, false);
    }


}


void compute_sigma_cdot_spinor( std::vector<mrcpp::CompFunction<3>*> &Curl_top, std::vector<mrcpp::CompFunction<3>*> &Curl_bottom, CompFunction<3> &SOC_Psi_t, CompFunction<3> &SOC_Psi_b){
    // Define the sigma matrices in a vector:
    std::vector<Eigen::Matrix2cd> sigma(3);
    sigma[0] << 0, 1,
                1, 0; // Pauli matrix sigma_x
    sigma[1] << 0, -std::complex<double>(0, 1),
                std::complex<double>(0, 1), 0; // Pauli matrix sigma_y
    sigma[2] << 1, 0,
                0, -1; // Pauli matrix sigma_z


    // I'll devide the xyz components of the top and bottoom components of the spinor in 3 spinors, one for each direction
    std::vector<mrcpp::CompFunction<3>> Spinor_component_x(2);
    std::vector<mrcpp::CompFunction<3>> Spinor_component_y(2);
    std::vector<mrcpp::CompFunction<3>> Spinor_component_z(2);

    Spinor_component_x[0] = *Curl_top[0];
    Spinor_component_x[1] = *Curl_bottom[0];

    Spinor_component_y[0] = *Curl_top[1];
    Spinor_component_y[1] = *Curl_bottom[1];

    Spinor_component_z[0] = *Curl_top[2];
    Spinor_component_z[1] = *Curl_bottom[2];


    // Now i compute the scalar product between the rotor and the sigma matrix vector
    // I'll do one by one
    // Compute SOC_i as the matrix-vector product of sigma_i and Spinor_component_i
    std::vector<mrcpp::CompFunction<3>*> SOC_x(2);
    std::vector<mrcpp::CompFunction<3>*> SOC_y(2);
    std::vector<mrcpp::CompFunction<3>*> SOC_z(2);

    // In general, for each matrix vector(spinor) multiplication we need 4 scalar multiplication and 2 additions

    // Sigma_x * Spinor_component_x
    mrcpp::add(*SOC_x[0], sigma[0](0,0), Spinor_component_x[0], sigma[0](0,1), Spinor_component_x[1], building_precision, false);
    mrcpp::add(*SOC_x[1], sigma[0](1,0), Spinor_component_x[0], sigma[0](1,1), Spinor_component_x[1], building_precision, false);

    // Sigma_y * Spinor_component_y
    mrcpp::add(*SOC_y[0], sigma[1](0,0), Spinor_component_y[0], sigma[1](0,1), Spinor_component_y[1], building_precision, false);
    mrcpp::add(*SOC_y[1], sigma[1](1,0), Spinor_component_y[0], sigma[1](1,1), Spinor_component_y[1], building_precision, false);

    // Sigma_z * Spinor_component_z
    mrcpp::add(*SOC_z[0], sigma[2](0,0), Spinor_component_z[0], sigma[2](0,1), Spinor_component_z[1], building_precision, false);
    mrcpp::add(*SOC_z[1], sigma[2](1,0), Spinor_component_z[0], sigma[2](1,1), Spinor_component_z[1], building_precision, false);


    // Now i sum the y and z components
    std::vector<mrcpp::CompFunction<3>> SOC_yz(2); // This is the sum of the y and z components, temporary, to be used in the final sum
    mrcpp::add(SOC_yz[0], 1.0, *SOC_y[0], 1.0, *SOC_z[0],building_precision, false);
    mrcpp::add(SOC_yz[1], 1.0, *SOC_y[1], 1.0, *SOC_z[1],building_precision, false);

    // Finally i sum the x component to the yz component to the total SOC
    mrcpp::add(SOC_Psi_t, 1.0, *SOC_x[0], 1.0, SOC_yz[0],building_precision, false);
    mrcpp::add(SOC_Psi_b, 1.0, *SOC_x[1], 1.0, SOC_yz[1],building_precision, false);

}





ComplexDouble compute_Term3_T_ZORA(MultiResolutionAnalysis<3> &MRA, std::vector<std::vector<mrcpp::CompFunction<3>*>> &Nabla_Psi_2c, mrcpp::CompFunction<3> &K_tree, std::vector<mrcpp::CompFunction<3> *> &Nabla_K_tree,  std::vector<mrcpp::CompFunction<3>> Psi_2c){
    // Compute the product K * (\Nabla \Psi) for top and bottom components
    std::vector<mrcpp::CompFunction<3>*> K_Nabla_Psi_top;
    std::vector<mrcpp::CompFunction<3>*> K_Nabla_Psi_bottom;

    // For synthax:
    // void mrcpp::multiply<3>(mrcpp::CompFunction<3> &out, mrcpp::CompFunction<3> inp_a, mrcpp::CompFunction<3> inp_b, double prec, bool absPrec, bool useMaxNorms, bool conjugate)
    for (int i = 0; i<2; i++){
        mrcpp::multiply(*K_Nabla_Psi_top[i], K_tree, *Nabla_Psi_2c[0][i] ,building_precision, false, false, false);   
        mrcpp::multiply(*K_Nabla_Psi_top[i], K_tree, *Nabla_Psi_2c[1][i] ,building_precision, false, false, false);   
    }
    
    // Compute the rotor of K * (\Nabla \Psi) for top and bottom components
    std::vector<mrcpp::CompFunction<3>*> rotor_K_Nabla_Psi_top;
    std::vector<mrcpp::CompFunction<3>*> rotor_K_Nabla_Psi_bottom;

    mrcpp::ABGVOperator<3> D(MRA, 0.0, 0.0);
    compute_rotor(rotor_K_Nabla_Psi_top, D, K_Nabla_Psi_top, MRA);
    compute_rotor(rotor_K_Nabla_Psi_bottom, D, K_Nabla_Psi_bottom, MRA);

    // Compute the scalar product between the rotor and the sigma matrix vector
    CompFunction<3> SOC_Psi_t(MRA);
    CompFunction<3> SOC_Psi_b(MRA);
    compute_sigma_cdot_spinor(rotor_K_Nabla_Psi_top, rotor_K_Nabla_Psi_bottom, SOC_Psi_t, SOC_Psi_b);

    // Now i do the innetr produc with the bra
    ComplexDouble Top_contribution = dot(Psi_2c[0], SOC_Psi_t);
    ComplexDouble Bottom_contribution = dot(Psi_2c[1], SOC_Psi_b);

    return ComplexDouble(0,1)*(Top_contribution + Bottom_contribution);
}







ComplexDouble compute_energy_ZORA(MultiResolutionAnalysis<3> &MRA, std::vector<std::vector<mrcpp::CompFunction<3>*>> &Nabla_Psi_2c, mrcpp::CompFunction<3> &K_tree, std::vector<mrcpp::CompFunction<3> *> &Nabla_K_tree,  std::vector<mrcpp::CompFunction<3>> Psi_2c, CompFunction<3> &V){
    ComplexDouble Term1 = compute_Term1_T_ZORA(MRA, Nabla_Psi_2c, K_tree, Psi_2c);
    ComplexDouble Term2 = compute_Term2_T_ZORA(MRA, Nabla_Psi_2c, K_tree, Nabla_K_tree, Psi_2c);
    ComplexDouble Term3 = compute_Term3_T_ZORA(MRA, Nabla_Psi_2c, K_tree, Nabla_K_tree, Psi_2c);

    // Now we do <Psi | V | Psi>
    mrcpp::CompFunction<3> psi_top__V(MRA);
    mrcpp::CompFunction<3> psi_bottom__V(MRA);

    mrcpp::multiply(psi_top__V, V, Psi_2c[0],building_precision, false, false, false);
    mrcpp::multiply(psi_bottom__V, V, Psi_2c[1],building_precision, false, false, false);

    ComplexDouble potential_energy = dot(Psi_2c[0],psi_top__V) + dot(Psi_2c[1],psi_bottom__V);

    ComplexDouble kinetic_energy = -0.5*m* (Term1 + Term2 + Term3) ;

    return kinetic_energy + potential_energy;



}
    









































double energy_ZORA(MultiResolutionAnalysis<3> &MRA, std::vector<mrcpp::CompFunction<3>> &Psi, CompFunction<3> &V){
    // Declare the FunctionTrees
    CompFunction<3> psi_top__V(MRA);
    CompFunction<3> psi_bottom__V(MRA);

    CompFunction<3> K_tree(MRA);
    CompFunction<3> K_psi_top_tree(MRA);

    // Define the K function
    std::function<double(const Coord<3> &x)> K_r = [] (const mrcpp::Coord<3> &r) -> double {
        double abs_r = std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        double result = 1.0 - 1.0/ ( abs_r * 2.0 * m * c * c);
        return 1.0 / result;
    };

    // And project it on the correspondig tree
    mrcpp::project(K_tree, K_r,building_precision);



    // Define the type of derivative operator: ABGV
    mrcpp::ABGVOperator<3> D(MRA, 0.0, 0.0);


    // *********************************************
    // Now we can start the subroutine:
 
 

    // ================================================
    // ============ KINETIC PART ====================
    // ================================================

    // TERM 1: <\Psi (\Nabla K) | (\Nabla \Psi) >  , for both the top and bottom components
    
    std::vector<mrcpp::CompFunction<3>> Psi_t_Nabla_K; 
    std::vector<mrcpp::CompFunction<3>> Psi_b_Nabla_K;

    CompFunction<3> Psi_Nabla_k_tmp(MRA);
    CompFunction<3> Psi_element(MRA);
    CompFunction<3> Nabla_K_element(MRA);
    
    

    // Gradient of K
    auto Nabla_K = mrcpp::gradient(D,K_tree);

    // Gradient of the 2 components of Psi
    // PAY ATTENTION!!! The gradient obtained is actually a pointer, so later on we'll have to dereference it
    auto Nabla_Psi_t = mrcpp::gradient(D,Psi[0]);
    auto Nabla_Psi_b = mrcpp::gradient(D,Psi[1]);

    // Compute Psi_top * Nabla_K
    for (int i=0;i<3;i++){
        mrcpp::multiply(Psi_Nabla_k_tmp, Psi[0], *Nabla_K[i],building_precision, false, false, true);   // <---- Here we dereference the pointer
        Psi_t_Nabla_K.push_back(Psi_Nabla_k_tmp);
    }
    
    // Compute Psi_bottom * Nabla_K
    for (int i=0;i<3;i++){
        mrcpp::multiply(Psi_Nabla_k_tmp, Psi[1], *Nabla_K[i],building_precision, false, false, true);   // <---- Here we dereference the pointer
        Psi_b_Nabla_K.push_back(Psi_Nabla_k_tmp);
    }

    std::vector<ComplexDouble> Kinetic_term_1;


    Kinetic_term_1.push_back(dot(Psi_t_Nabla_K[0],*Nabla_Psi_t[0]) + dot(Psi_t_Nabla_K[1],*Nabla_Psi_t[1]) + dot(Psi_t_Nabla_K[2],*Nabla_Psi_t[2]));
    Kinetic_term_1.push_back(dot(Psi_b_Nabla_K[0],*Nabla_Psi_b[0]) + dot(Psi_b_Nabla_K[1],*Nabla_Psi_b[1]) + dot(Psi_b_Nabla_K[2],*Nabla_Psi_b[2]));

    // TERM 2: - <\Nabla (\Psi K) | \Nabla \Psi>

    CompFunction<3> Psi_t_K;
    CompFunction<3> Psi_b_K;


    // Compute Psi_top * K
    mrcpp::multiply(Psi_t_K, Psi[0], K_tree,building_precision, false, false, true);
    // Compute Psi_bottom * K
    mrcpp::multiply(Psi_b_K, Psi[1], K_tree,building_precision, false, false, true);

    // Gradient of Psi_top * K
    auto Nabla_Psi_t_K = mrcpp::gradient(D,Psi_t_K);
    // Gradient of Psi_bottom * K
    auto Nabla_Psi_b_K = mrcpp::gradient(D,Psi_b_K);

    // Now we compute the contribution from term 2:
    std::vector<ComplexDouble> Kinetic_term_2;

    // REMEMBER THE MINUS SIGN!!! (it comes from the integration by parts)
    Kinetic_term_2.push_back(-1.0*(dot(*Nabla_Psi_t_K[0],*Nabla_Psi_t[0]) + dot(*Nabla_Psi_t_K[1],*Nabla_Psi_t[1]) + dot(*Nabla_Psi_t_K[2],*Nabla_Psi_t[2])));
    Kinetic_term_2.push_back(-1.0*(dot(*Nabla_Psi_b_K[0],*Nabla_Psi_b[0]) + dot(*Nabla_Psi_b_K[1],*Nabla_Psi_b[1]) + dot(*Nabla_Psi_b_K[2],*Nabla_Psi_b[2])));



    // TERM 3: < \Psi | i * \sigma \cdot (\Nabla \cross K (\Nabla)) | \Psi >

    // Define the sigma matrices in a vector:
        std::vector<Eigen::Matrix2cd> sigma(3);
        sigma[0] << 0, 1,
                    1, 0; // Pauli matrix sigma_x
        sigma[1] << 0, -std::complex<double>(0, 1),
                    std::complex<double>(0, 1), 0; // Pauli matrix sigma_y
        sigma[2] << 1, 0,
                    0, -1; // Pauli matrix sigma_z

    // Define the Levi-Civita tensor for 3D
    /* ComplexDouble levi_civita[3][3][3] = {
        {{0.0, 0.0, 0.0},
         {0.0, 0.0, 1.0},
         {0.0, -1.0, 0.0}},
        {{0.0, 0.0, -1.0}, 
        {0.0, 0.0, 0.0}, 
        {1.0, 0.0, 0.0}},
        {{0.0, 1.0, 0.0}, 
        {-1.0, 0.0, 0.0}, 
        {0.0, 0.0, 0.0}}
    };
 */
    ComplexDouble levi_civita[3][3][3] = {
        {{ComplexDouble(0.0, 0.0), ComplexDouble(0.0, 0.0), ComplexDouble(0.0, 0.0)}, 
         {ComplexDouble(0.0, 0.0), ComplexDouble(0.0, 0.0), ComplexDouble(1.0, 0.0)}, 
         {ComplexDouble(0.0, 0.0), ComplexDouble(-1.0, 0.0), ComplexDouble(0.0, 0.0)}},
        {{ComplexDouble(0.0, 0.0), ComplexDouble(0.0, 0.0), ComplexDouble(-1.0, 0.0)}, 
         {ComplexDouble(0.0, 0.0), ComplexDouble(0.0, 0.0), ComplexDouble(0.0, 0.0)}, 
         {ComplexDouble(1.0, 0.0), ComplexDouble(0.0, 0.0), ComplexDouble(0.0, 0.0)}},
        {{ComplexDouble(0.0, 0.0), ComplexDouble(1.0, 0.0), ComplexDouble(0.0, 0.0)}, 
         {ComplexDouble(-1.0, 0.0), ComplexDouble(0.0, 0.0), ComplexDouble(0.0, 0.0)}, 
         {ComplexDouble(0.0, 0.0), ComplexDouble(0.0, 0.0), ComplexDouble(0.0, 0.0)}}
    };
   
    // Compute the cross product of Nabla and K
    std::vector<CompFunction<3>> SOC_Psi(2);


    CompFunction<3> Nabla_cross_K__Nabla_Psi_tmp(MRA);
    //std::vector<CompFunction<3>> deriv_Psi;
    
    CompFunction<3> deriv_Psi1(MRA);
    CompFunction<3> deriv_K1(MRA);
    CompFunction<3> deriv_Psi2(MRA);
    CompFunction<3> deriv_K2(MRA);

    CompFunction<3> Rotor_component1(MRA);
    CompFunction<3> Rotor_component2(MRA);

    std::vector<std::vector<CompFunction<3>>> rotor_K_Nabla_Psi(2, std::vector<CompFunction<3>>(3, CompFunction<3>(MRA)));

    for (int s = 0; s<2;s++){
        for (int k=0; k<3; k++){
            for (int i=0; i<3; i++){
                // Skip if i == k because the Levicivita tensor is antisymmetric
                if (i == k){
                    continue;
                }
                
                
                for (int j=i+1; j<3;j++){
                // Skip if j == i or j == k because the Levicivita tensor is antisymmetric
                    if (j==i || j==k){
                        continue;
                    }

                    mrcpp::apply(deriv_Psi1, D, Psi[s], j); // Derivative of Psi with respect to j
                    mrcpp::apply(deriv_K1, D, K_tree, i);   // Derivative of K with respect to i
                    mrcpp::apply(deriv_Psi2, D, Psi[s], i); // Derivative of Psi with respect to i
                    mrcpp::apply(deriv_K2, D, K_tree, j);   // Derivative of K with respect to j

                    // I create the 2 terms of the rotor, that is the moltiplication between the derivatives for both Psi and K in the 2 remaining directions
                    mrcpp::multiply(Rotor_component1, deriv_Psi1, deriv_K1,building_precision, false, false, false);
                    mrcpp::multiply(Rotor_component2, deriv_Psi2, deriv_K2,building_precision, false, false, false);

                    // I compute the rotor by adding the 2 terms with the correct sign given by the Levi-Civita tensor
                    mrcpp::add<3>(rotor_K_Nabla_Psi[s][k], levi_civita[k][i][j], Rotor_component1, levi_civita[k][j][i], Rotor_component2, building_precision, false);
                }
            }        
        }
    }

    // Now in rotor_K_Nabla_Psi i have the 3 components of the rotor (k = [0,2]) for both the components of the spinor (s = [0,1])
    // Let's compute the scalar product between the rotor and the sigma matrix vector
    // In order not to get confused i'll compute the scalar product between the rotor and the sigma matrix vector in a separate loop

    /*
     * rotor_K_Nabla_Psi = { Rot(K*Grad*PsiTop)_x, Rot(K*Grad*PsiTop)_y, Rot(K*Grad*PsiTop)_z; 
     *                       Rot(K*Grad*PsiBottom)_x, Rot(K*Grad*PsiBottom)_y, Rot(K*Grad*PsiBottom)_z}  
    */

    std::vector<CompFunction<3>> SOC_Psi_x(2);
    std::vector<CompFunction<3>> SOC_Psi_y(2);
    std::vector<CompFunction<3>> SOC_Psi_z(2);

    std::vector<CompFunction<3>> SOC_Psi_yz(2); // This is the sum of the y and z components, temporary, to be used in the final sum


    for (int r = 0; r<2; r++){
        // Sigma_z * rotor_K_Nabla_Psi_z    
        mrcpp::add(SOC_Psi_z[r], sigma[2](r,0), rotor_K_Nabla_Psi[0][2], sigma[2](r,1), rotor_K_Nabla_Psi[1][2], building_precision, false);

        // Sigma_y * rotor_K_Nabla_Psi_y
        mrcpp::add(SOC_Psi_y[r], sigma[1](r,0), rotor_K_Nabla_Psi[0][1], sigma[1](r,1), rotor_K_Nabla_Psi[1][1], building_precision, false);

        // Sigma_x * rotor_K_Nabla_Psi_x
        mrcpp::add(SOC_Psi_x[r], sigma[0](r,0), rotor_K_Nabla_Psi[0][0], sigma[0](r,1), rotor_K_Nabla_Psi[1][0], building_precision, false);


        // Now i sum the y and z components
        mrcpp::add(SOC_Psi_yz[r], 1.0, SOC_Psi_y[r], 1.0, SOC_Psi_z[r],building_precision, false);

        // Finally i sum the x component to the yz component to the total SOC
        mrcpp::add(SOC_Psi[r], 1.0, SOC_Psi_x[r], 1.0, SOC_Psi_yz[r],building_precision, false);               
    }
    
    std::vector<ComplexDouble> Kinetic_term_3;
    SOC_Psi[0].rescale(ComplexDouble(0.0, 1.0));
    SOC_Psi[1].rescale(ComplexDouble(0.0, 1.0));

    Kinetic_term_3.push_back(dot( Psi[0], SOC_Psi[0]));
    Kinetic_term_3.push_back(dot( Psi[1], SOC_Psi[1]));



    // ================================================
    // ============ POTENTIAL PART ====================
    // ================================================
    mrcpp::multiply(building_precision, psi_top__V, 1.0, Psi[0], V);
    mrcpp::multiply(building_precision, psi_bottom__V, 1.0, Psi[1], V);
    ComplexDouble potential_energy = dot(Psi[0],psi_top__V) + dot(Psi[1],psi_bottom__V);


    // ================================================
    // ============ TOTAL ENERGY ====================
    // ================================================

    ComplexDouble tot_kin = Kinetic_term_1[0] + Kinetic_term_1[1] + Kinetic_term_2[0] + Kinetic_term_2[1] + Kinetic_term_3[0] + Kinetic_term_3[1];
    tot_kin = tot_kin / (2.0*m);

    ComplexDouble result = potential_energy + tot_kin;

    std::cout << "Real energy = " << result.real() << '\n';
    std::cout << "Imag energy = " << result.imag() << '\n';


    return result.real();
}



void apply_Helmholtz_ZORA(MultiResolutionAnalysis<3> &MRA, std::vector<mrcpp::CompFunction<3>> &Psi_in, std::vector<mrcpp::CompFunction<3>> &Psi_out, CompFunction<3> &V, double mu){
    // Build the Helmholtz operator
    mrcpp::HelmholtzOperator Helm(MRA, mu, building_precision);
    double energy = - (mu*mu)/(2*m);

    // we will apply the Helmholtz operator to the 2 components of the spinor separately
    // Psi_out = Helm(-term_1 + K_inverted * (Kinetic_term_2 + Kinetic_term_3 + 2*m*V)  )

    /*
     * Where the term_1 = (energy/c^2) * V * Psi_in
     * K_inverted = (1 - V / (2*m*c^2) = 1/K_r
     * Kinetic_term_2 =  \Nabla(K) \cdot \Nabla(\Psi_in)
     * Kinetic_term_3 = i * \sigma \cdot (\Nabla \cross K (\Nabla (\Psi_in)))
     *     
     * 
    */


    // ================================================
    // ================= TERM 1 =======================
    // ================================================
    // term_1 = (energy/c^2) * V * Psi_in

    CompFunction<3> term_1_top(MRA);
    CompFunction<3> term_1_bottom(MRA);

    mrcpp::multiply(term_1_top, V, Psi_in[0],building_precision, false, false, false);
    mrcpp::multiply(term_1_bottom, V, Psi_in[1],building_precision, false, false, false);

    double scaling_factor_term_1 = energy / (c*c);

    term_1_top.rescale(scaling_factor_term_1);
    term_1_bottom.rescale(scaling_factor_term_1);


    // ================================================
    // ================= K_inverted ===================
    // ================================================

    // Define the K function
    std::function<double(const Coord<3> &x)> K_r = [] (const mrcpp::Coord<3> &r) -> double {
        double abs_r = std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        double result = 1.0 - (1.0/ ( abs_r * 2.0 * m * c * c));
        return 1.0 / result;
    };

    std::function<double(const Coord<3> &x)> K_r_inverted = [] (const mrcpp::Coord<3> &r) -> double {
        double abs_r = std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        double result = 1.0 - (1.0/ ( abs_r * 2.0 * m * c * c));
        return result;
    };

    // Define the type of derivative operator: ABGV
    mrcpp::ABGVOperator<3> D(MRA, 0.0, 0.0);

    // Define 1/K function
    CompFunction<3> K_inverted(MRA);
    CompFunction<3> K_tree(MRA);
    mrcpp::project(K_inverted, K_r_inverted,building_precision);
    mrcpp::project(K_tree, K_r,building_precision);

   
    // ================================================
    // ================ 2 * m * V =====================
    // ================================================
    // two_m_V = 2 * m * V * Psi_in

    CompFunction<3> two_m_V_top(MRA);
    CompFunction<3> two_m_V_bottom(MRA);

    mrcpp::multiply(two_m_V_top,      V, Psi_in[0], building_precision, false, false, false);
    mrcpp::multiply(two_m_V_bottom,   V, Psi_in[1], building_precision, false, false, false);

    two_m_V_top.rescale(ComplexDouble(2.0*m,0.0));
    two_m_V_bottom.rescale(ComplexDouble(2.0*m,0.0));



    // ================================================
    // ============= Kinetic term 2 ===================
    // ================================================

    // Kinetic_term_2 =  \Nabla(K) \cdot \Nabla(\Psi_in)

    // Define the gradient of K
    auto Nabla_K = mrcpp::gradient(D,K_tree);


    // Gradient of the 2 components of Psi
    // PAY ATTENTION!!! The gradient obtained is actually a pointer, so later on we'll have to dereference it
    auto Nabla_Psi_t = mrcpp::gradient(D,Psi_in[0]);
    auto Nabla_Psi_b = mrcpp::gradient(D,Psi_in[1]);


    // Compute Nabla_K * Nabla_Psi for each component
    
    // z component:
    mrcpp::CompFunction<3> Kinetic_term_2_top_z(MRA);
    mrcpp::CompFunction<3> Kinetic_term_2_bottom_z(MRA);

    mrcpp::multiply(Kinetic_term_2_top_z, *Nabla_K[2], *Nabla_Psi_t[2],building_precision, false, false, false);
    mrcpp::multiply(Kinetic_term_2_bottom_z, *Nabla_K[2], *Nabla_Psi_b[2],building_precision, false, false, false);
    
    // y component:
    mrcpp::CompFunction<3> Kinetic_term_2_top_y(MRA);
    mrcpp::CompFunction<3> Kinetic_term_2_bottom_y(MRA);

    mrcpp::multiply(Kinetic_term_2_top_y, *Nabla_K[1], *Nabla_Psi_t[1],building_precision, false, false, false);
    mrcpp::multiply(Kinetic_term_2_bottom_y, *Nabla_K[1], *Nabla_Psi_b[1],building_precision, false, false, false);

    // x component:
    mrcpp::CompFunction<3> Kinetic_term_2_top_x(MRA);
    mrcpp::CompFunction<3> Kinetic_term_2_bottom_x(MRA);

    mrcpp::multiply(Kinetic_term_2_top_x, *Nabla_K[0], *Nabla_Psi_t[0],building_precision, false, false, false);
    mrcpp::multiply(Kinetic_term_2_bottom_x, *Nabla_K[0], *Nabla_Psi_b[0],building_precision, false, false, false);

    // Now i sum the 3 components
    // These are the first two additions
    mrcpp::CompFunction<3> Kinetic_term_2_top_zy(MRA);
    mrcpp::CompFunction<3> Kinetic_term_2_bottom_zy(MRA);

    // Then the total, by adding x as well
    mrcpp::CompFunction<3> Kinetic_term_2_top(MRA);
    mrcpp::CompFunction<3> Kinetic_term_2_bottom(MRA);
       
    mrcpp::add(Kinetic_term_2_top_zy,    1.0,    Kinetic_term_2_top_z,    1.0,       Kinetic_term_2_top_y,       building_precision, false);
    mrcpp::add(Kinetic_term_2_bottom_zy, 1.0,    Kinetic_term_2_bottom_z, 1.0,       Kinetic_term_2_bottom_y,    building_precision, false);
    
    mrcpp::add(Kinetic_term_2_top,       1.0,    Kinetic_term_2_top_zy,    1.0,       Kinetic_term_2_top_x,       building_precision, false);
    mrcpp::add(Kinetic_term_2_bottom,    1.0,    Kinetic_term_2_bottom_zy, 1.0,    Kinetic_term_2_bottom_x,       building_precision, false);
    
    // ================================================
    // ============= Kinetic term 3 ===================
    // ================================================
    // Kinetic_term_3 = i * \sigma \cdot (\Nabla \cross K (\Nabla (\Psi_in)))

    // Define the sigma matrices in a vector:
    std::vector<Eigen::Matrix2cd> sigma(3);
    sigma[0] << 0, 1,
                1, 0; // Pauli matrix sigma_x
    sigma[1] << 0, -std::complex<double>(0, 1),
                std::complex<double>(0, 1), 0; // Pauli matrix sigma_y
    sigma[2] << 1, 0,
                0, -1; // Pauli matrix sigma_z      

    // Define the Levi-Civita tensor for 3D
    ComplexDouble levi_civita[3][3][3] = {
        {{ComplexDouble(0.0, 0.0), ComplexDouble(0.0, 0.0), ComplexDouble(0.0, 0.0)}, 
         {ComplexDouble(0.0, 0.0), ComplexDouble(0.0, 0.0), ComplexDouble(1.0, 0.0)}, 
         {ComplexDouble(0.0, 0.0), ComplexDouble(-1.0, 0.0), ComplexDouble(0.0, 0.0)}},
        {{ComplexDouble(0.0, 0.0), ComplexDouble(0.0, 0.0), ComplexDouble(-1.0, 0.0)}, 
         {ComplexDouble(0.0, 0.0), ComplexDouble(0.0, 0.0), ComplexDouble(0.0, 0.0)}, 
         {ComplexDouble(1.0, 0.0), ComplexDouble(0.0, 0.0), ComplexDouble(0.0, 0.0)}},
        {{ComplexDouble(0.0, 0.0), ComplexDouble(1.0, 0.0), ComplexDouble(0.0, 0.0)}, 
         {ComplexDouble(-1.0, 0.0), ComplexDouble(0.0, 0.0), ComplexDouble(0.0, 0.0)}, 
         {ComplexDouble(0.0, 0.0), ComplexDouble(0.0, 0.0), ComplexDouble(0.0, 0.0)}}
    };
   
    // Compute the cross product of Nabla and K
    std::vector<CompFunction<3>> SOC_Psi(2);


    CompFunction<3> Nabla_cross_K__Nabla_Psi_tmp(MRA);
    
    
    CompFunction<3> deriv_Psi1(MRA);
    CompFunction<3> deriv_K1(MRA);
    CompFunction<3> deriv_Psi2(MRA);
    CompFunction<3> deriv_K2(MRA);

    CompFunction<3> Rotor_component1(MRA);
    CompFunction<3> Rotor_component2(MRA);

    std::vector<std::vector<CompFunction<3>>> rotor_K_Nabla_Psi(2, std::vector<CompFunction<3>>(3, CompFunction<3>(MRA)));

    for (int s = 0; s<2;s++){
        for (int k=0; k<3; k++){
            for (int i=0; i<3; i++){
                // Skip if i == k because the Levicivita tensor is antisymmetric
                if (i == k){
                    continue;
                }
                
                
                for (int j=i+1; j<3;j++){
                // Skip if j == i or j == k because the Levicivita tensor is antisymmetric
                    if (j==i || j==k){
                        continue;
                    }

                    mrcpp::apply(deriv_Psi1, D, Psi_in[s], j); // Derivative of Psi with respect to j
                    mrcpp::apply(deriv_K1, D, K_tree, i);   // Derivative of K with respect to i
                    mrcpp::apply(deriv_Psi2, D, Psi_in[s], i); // Derivative of Psi with respect to i
                    mrcpp::apply(deriv_K2, D, K_tree, j);   // Derivative of K with respect to j

                    // I create the 2 terms of the rotor, that is the moltiplication between the derivatives for both Psi and K in the 2 remaining directions
                    mrcpp::multiply(Rotor_component1, deriv_Psi1, deriv_K1,building_precision, false, false, false);
                    mrcpp::multiply(Rotor_component2, deriv_Psi2, deriv_K2,building_precision, false, false, false);

                    // I compute the rotor by adding the 2 terms with the correct sign given by the Levi-Civita tensor
                    mrcpp::add<3>(rotor_K_Nabla_Psi[s][k], levi_civita[k][i][j], Rotor_component1, levi_civita[k][j][i], Rotor_component2, building_precision, false);
                }
            }        
        }
    }

    // Now in rotor_K_Nabla_Psi i have the 3 components of the rotor (k = [0,2]) for both the components of the spinor (s = [0,1])
    // Let's compute the scalar product between the rotor and the sigma matrix vector
    // In order not to get confused i'll compute the scalar product between the rotor and the sigma matrix vector in a separate loop

    /*
     * rotor_K_Nabla_Psi = { Rot(K*Grad*PsiTop)_x, Rot(K*Grad*PsiTop)_y, Rot(K*Grad*PsiTop)_z; 
     *                       Rot(K*Grad*PsiBottom)_x, Rot(K*Grad*PsiBottom)_y, Rot(K*Grad*PsiBottom)_z}  
    */

    std::vector<CompFunction<3>> SOC_Psi_x(2);
    std::vector<CompFunction<3>> SOC_Psi_y(2);
    std::vector<CompFunction<3>> SOC_Psi_z(2);

    std::vector<CompFunction<3>> SOC_Psi_yz(2); // This is the sum of the y and z components, temporary, to be used in the final sum


    for (int r = 0; r<2; r++){
        // Sigma_z * rotor_K_Nabla_Psi_z    
        mrcpp::add<3>(SOC_Psi_z[r], sigma[2](r,0), rotor_K_Nabla_Psi[0][2], sigma[2](r,1), rotor_K_Nabla_Psi[1][2], building_precision, false);

        // Sigma_y * rotor_K_Nabla_Psi_y
        mrcpp::add<3>(SOC_Psi_y[r], sigma[1](r,0), rotor_K_Nabla_Psi[0][1], sigma[1](r,1), rotor_K_Nabla_Psi[1][1], building_precision, false);

        // Sigma_x * rotor_K_Nabla_Psi_x
        mrcpp::add<3>(SOC_Psi_x[r], sigma[0](r,0), rotor_K_Nabla_Psi[0][0], sigma[0](r,1), rotor_K_Nabla_Psi[1][0], building_precision, false);


        // Now i sum the y and z components
        mrcpp::add<3>(SOC_Psi_yz[r], 1.0, SOC_Psi_y[r], 1.0, SOC_Psi_z[r],building_precision, false);

        // Finally i sum the x component to the yz component to the total SOC
        mrcpp::add<3>(SOC_Psi[r], 1.0, SOC_Psi_x[r], 1.0, SOC_Psi_yz[r],building_precision, false);               
    }

    SOC_Psi[0].rescale(ComplexDouble(0.0, 1.0));
    SOC_Psi[1].rescale(ComplexDouble(0.0, 1.0));

    // ================================================
    // =============== TOTAL TERM 2 ===================
    // ================================================
    // Total_term_2 = K_inverted * (Kinetic_term_2 + Kinetic_term_3 - 2*m*V) 


    CompFunction<3> Total_term_2_top(MRA);
    CompFunction<3> Total_term_2_bottom(MRA);

    CompFunction<3> Total_term_2_top_addition_tmp(MRA);
    CompFunction<3> Total_term_2_bottom_addition_tmp(MRA);

    CompFunction<3> Total_term_2_top_addition(MRA);
    CompFunction<3> Total_term_2_bottom_addition(MRA);


    // first two additions
    mrcpp::add<3>(Total_term_2_top_addition_tmp, 1.0, Kinetic_term_2_top, 1.0, SOC_Psi[0], building_precision, false);
    mrcpp::add<3>(Total_term_2_bottom_addition_tmp, 1.0, Kinetic_term_2_bottom, 1.0, SOC_Psi[1], building_precision, false);

    // then the total, by adding 2*m*V as well
    mrcpp::add<3>(Total_term_2_top_addition, 1.0,   Total_term_2_top_addition_tmp, -1.0, two_m_V_top, building_precision, false);
    mrcpp::add<3>(Total_term_2_bottom_addition,  1.0,   Total_term_2_bottom_addition_tmp, -1.0, two_m_V_bottom, building_precision, false);

    // Now i multiply by K_inverted
    mrcpp::multiply(Total_term_2_top, K_inverted, Total_term_2_top_addition, building_precision, false, false, false);
    mrcpp::multiply(Total_term_2_bottom, K_inverted, Total_term_2_bottom_addition, building_precision, false, false, false);
   
    // ================================================
    // =================== Final ======================
    // ================================================

    // V_ZORA_Psi = -term_1 + Total_term_2

    CompFunction<3> V_ZORA_Psi_top(MRA);
    CompFunction<3> V_ZORA_Psi_bottom(MRA);

    mrcpp::add<3>(V_ZORA_Psi_top, -1.0, term_1_top, 1.0, Total_term_2_top, building_precision, false);
    mrcpp::add<3>(V_ZORA_Psi_bottom, -1.0, term_1_bottom, 1.0, Total_term_2_bottom, building_precision, false);



    // Apply the Helmholtz operator
    mrcpp::apply(building_precision,Psi_out[0], Helm, Psi_in[0]);
    mrcpp::apply(building_precision,Psi_out[1], Helm, Psi_in[1]);

}



