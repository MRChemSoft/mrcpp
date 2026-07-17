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


// THIS IS ONLY FOR HELIUM

void Update_SCF_Variables(MultiResolutionAnalysis<3> &MRA, ABGVOperator<3> &D, std::vector<mrcpp::CompFunction<3>> &Psi_2c, CompFunction<3> &V_electron_nucleus ,std::vector<std::vector<mrcpp::CompFunction<3> *>> &Nabla_Psi_2c,CompFunction<3> &K_tree, std::vector<mrcpp::CompFunction<3> *> &Nabla_K_tree, CompFunction<3> &K_inverse, CompFunction<3> &V, CompFunction<3> &J, mrcpp::PoissonOperator &P){
    // Updating the Nabla_Psi_2c
    //std::cout << "--------------UPDATING GRADIENT--------------" << '\n';
    Nabla_Psi_2c[0] = mrcpp::gradient(D,Psi_2c[0]);
    Nabla_Psi_2c[1] = mrcpp::gradient(D,Psi_2c[1]);
    // Print  values
    if (debug){
    std::cout << "Nabla_Psi_2c[0] = " << Nabla_Psi_2c[0][0]->getSquareNorm() << "   " <<  Nabla_Psi_2c[0][1]->getSquareNorm() << "   " <<  Nabla_Psi_2c[0][2]->getSquareNorm() << '\n';
    std::cout << "Nabla_Psi_2c[1] = " << Nabla_Psi_2c[1][0]->getSquareNorm() << "   " <<  Nabla_Psi_2c[1][1]->getSquareNorm() << "   " <<  Nabla_Psi_2c[1][2]->getSquareNorm() << '\n';
    std::cout << '\n';
    
    // Same format, print the pointers
    std::cout << "Nabla_Psi_2c[0][0] = " << Nabla_Psi_2c[0][0] << "   " << Nabla_Psi_2c[0][1] << "   " << Nabla_Psi_2c[0][2] << '\n';
    std::cout << "Nabla_Psi_2c[1][0] = " << Nabla_Psi_2c[1][0] << "   " << Nabla_Psi_2c[1][1] << "   " << Nabla_Psi_2c[1][2] << '\n';
    }
    //std::cout << "------------------- DONE -------------------" << '\n' << '\n';
    
    // Create the Rho_t and Rho_b, then Rho = Rho_t + Rho_b
    //std::cout << "--------------UPDATING RHO------------------" << '\n';
    CompFunction<3> Rho_t(MRA);
    CompFunction<3> Rho_b(MRA);
    CompFunction<3> Rho(MRA);


    //std::cout << "Computing density for TOP component..." << '\n';
    make_density_local(Rho_t, Psi_2c[0] , MRA, building_precision);
    //std::cout << "Norm of Rho_top = " << Rho_t.getSquareNorm() << '\n';

  
    //std::cout << "Computing density for BOTTOM component..." << '\n';
    make_density_local(Rho_b, Psi_2c[1] , MRA, building_precision);
    //std::cout << "Norm of Rho_bottom = " << Rho_b.getSquareNorm() << '\n' << '\n';
  

    mrcpp::add(Rho, 1.0, Rho_t, 1.0, Rho_b, building_precision, false);
    
    if ( num_cycle !=0  && (Rho_b.getSquareNorm() == 0.0 || Rho_t.getSquareNorm() == 0.0)) {
        std::cerr << "Error: Rho_t or Rho_b has zero norm!" << '\n';
        exit(1);
        return;
    }

    //std::cout << "---------------- DONE -------------------" << '\n' << '\n';

    // Calculate the J operator
    //std::cout << "--------------UPDATING J-----------------" << '\n';
    //std::cout << "Rho status" << *(Rho.CompD[0]) << '\n';
    mrcpp::apply(building_precision, J, P, Rho);
    //std::cout << "Applied the Poisson operator to the density" << '\n';
    double two_pi = (2.0 * M_PI);
    J.rescale(ComplexDouble(2.0 ,0.0));  
    //std::cout << "---------------- DONE -------------------" << '\n' << '\n';


    // Compute the full potential
    //std::cout << "--------------UPDATING V------------------" << '\n';
    mrcpp::add(V, 1.0, V_electron_nucleus,0.5, J, building_precision, false); // V = V_electron_nucleus + 0.5 J = V_ce(r)+J(r)
    //std::cout << "V norm = " << V.getSquareNorm() << '\n';
    //std::cout << "---------------- DONE -------------------" << '\n' << '\n';





    // Compute the inverse of K_inverse (that is just the Kappa operator):
    //std::cout << "--------------UPDATING K_tree------------------" << '\n';
    double constant = 2.0 * m * c * c;
    double prefactor = 1/constant;


    // Define and project a the function whcih is one everywhere
    std::function<double(const Coord<3> &x)> one = [] (const mrcpp::Coord<3> &r) -> double {
        return 1.0;
    };
    mrcpp::CompFunction<3> one_tree(MRA);
    mrcpp::project(one_tree, one, building_precision);

    // Function to calculate Kappa
    std::function<double(const Coord<3> &x)> K_function_compute = [&V_electron_nucleus, &J, prefactor, &V] (const mrcpp::Coord<3> &r) -> double {
        double VAL = 1.0 - (prefactor * V.CompD[0]->evalf(r));
        return (1.0 / VAL) -1.0;   
    };
    mrcpp::project(K_tree, K_function_compute, building_precision);
    //std::cout << "K_tree norm NO MAP= " << K_tree.getSquareNorm() << '\n';
 
    K_tree.add(ComplexDouble(1.0,0.0), one_tree); // K_tree = K_tree + 1

    
    if (debug){
    std::cout << "K_tree norm = " << K_tree.getSquareNorm() << '\n';
    std::cout << "K_tree is real = " << K_tree.isreal() << '\n';
    std::cout << "K_tree is complex = " << K_tree.iscomplex() << '\n';
    }
    //std::cout << "---------------- DONE -------------------" << '\n' << '\n';





    // Compute the new K tree
    //std::cout << "-------------UPDATING K_inv---------------" << '\n';


    // K_inverse
    std::function<double(const Coord<3> &x)> K_inverse_function = [&K_tree, &V, prefactor, &V_electron_nucleus, &J] (const mrcpp::Coord<3> &r) -> double {
        return - (prefactor * (V.CompD[0]->evalf(r)));
        //return 1.0 / K_tree.CompD[0]->evalf(r);   
    };

    mrcpp::project(K_inverse, K_inverse_function, building_precision);
    //add: K_inverse + 1
    K_inverse.add(ComplexDouble(1.0,0.0), one_tree);
    

    
    if (debug){
    std::cout << "K_inverse norm = " << K_inverse.getSquareNorm() << '\n';
    std::cout << "K_inverse is real = " << K_inverse.isreal() << '\n';
    std::cout << "K_inverse is complex = " << K_inverse.iscomplex() << '\n';
    std::cout << "K_inv Function Tree:" << '\n';
    std::cout << *(K_inverse.CompD[0]) << '\n';
    }
    



    //std::cout << "---------------- DONE -------------------" << '\n' << '\n';



    



    // Compute the Nabla_K_tree
    //std::cout << "--------------UPDATING Nabla_K_tree------------------" << '\n';
    Nabla_K_tree = mrcpp::gradient(D,K_tree);

    if (debug){
    // Print the pointers of the Nabla_K_tree
    std::cout << "Nabla_K_tree[0] = " << Nabla_K_tree[0] << '\n';
    std::cout << "Nabla_K_tree[1] = " << Nabla_K_tree[1] << '\n';
    std::cout << "Nabla_K_tree[2] = " << Nabla_K_tree[2] << '\n' <<'\n';
    // Now we print the norms
    std::cout << "Nabla_K_tree[0] norm = " << Nabla_K_tree[0]->getSquareNorm() << '\n';
    std::cout << "Nabla_K_tree[1] norm = " << Nabla_K_tree[1]->getSquareNorm() << '\n';
    std::cout << "Nabla_K_tree[2] norm = " << Nabla_K_tree[2]->getSquareNorm() << '\n' <<'\n';

    std::cout << "Nabla_K_tree[0] is real = " << Nabla_K_tree[0]->isreal() << '\n';
    std::cout << "Nabla_K_tree[0] is complex = " << Nabla_K_tree[0]->iscomplex() << '\n';
    }
    //std::cout << "-------------- DONE UPDATING SCF VARIABLES -------------------" << '\n' << '\n';

  
}