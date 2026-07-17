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







void make_density_local(CompFunction<3> &out, CompFunction<3> &inp, MultiResolutionAnalysis<3> &MRA, double prec) {
    //std::cout << '\t' << "make_density_local started, multiplying..." << '\n';

    CompFunction<3> inp_2(MRA);
    deep_copy(inp_2, inp);

    mrcpp::multiply(prec, out, 1.0, inp_2, inp, -1, false, false, true);
    if (debug){
    std::cout << '\t' << "out.iscomplex() = " << out.iscomplex() << '\n';
    std::cout << '\t' << "out.isreal() = " << out.isreal() << '\n'<<'\n';
    }
    
    if (out.iscomplex()) {
        FunctionTree<3, double> *out_tree_REAL = out.CompC[0]->Real();
        if (debug){
        std::cout << '\t' << "out.CompD[0] =" << out.CompD[0]<< '\n';
        std::cout << '\t' << "out.CompC[0] =" << out.CompC[0]<< '\n' << '\n';

        std::cout << '\t' << "*out.CompC[0].Real =" << out_tree_REAL << " <- We should be the same value!" << '\n';
        std::cout << '\t' << *out_tree_REAL << '\n';
 
        std::cout << '\t' << "out.CompC[0] =" << *(out.CompC[0])  << '\n';
        }
 
        out.CompD[0] = out_tree_REAL;
        out.CompC[0] = nullptr; // Set the pointer to nullptr to avoid double delete

        out.func_ptr->isreal = 1;
        out.func_ptr->iscomplex = 0;
    }
}



/*
*
*   CHANGE IN TERM 1 AND 2: THE FUNCTION DOT ALREADY CONJUGATES THE COMPFUNCTION BRA!!!!!!!!! UNCONJUGATE THEM AS IT IS DONE IN THE FUNCTION
*
*/


// IN principle, one could skip the step of recomputing all the gradient of the \Psi * K bu simpli adding \Nabla \Psi  K + \Psi  \Nabla K, but this would be longer, tho requiring less memory
ComplexDouble compute_Term1_T_ZORA(MultiResolutionAnalysis<3> &MRA, std::vector<std::vector<mrcpp::CompFunction<3>*>> &Nabla_Psi_2c, mrcpp::CompFunction<3> &K_tree,  std::vector<mrcpp::CompFunction<3>> &Psi_2c){
    // Nabla(\Psi K) = \Nabla \Psi  K + \Psi  \Nabla K
    debug = true;
    CompFunction<3> Psi_t_K(MRA);
    CompFunction<3> Psi_b_K(MRA);

    // Print the norm of the components it takes as input:
    if (debug){
    std::cout << "Psi_2c[0] = " << Psi_2c[0].getSquareNorm() << '\n';
    std::cout << "Psi_2c[1] = " << Psi_2c[1].getSquareNorm() << '\n';
    std::cout << "K_tree = " << K_tree.getSquareNorm() << '\n';
    std::cout << "Nabla_Psi_2c[0][0] = " << Nabla_Psi_2c[0][0]->getSquareNorm() << '\n';
    std::cout << "Nabla_Psi_2c[0][1] = " << Nabla_Psi_2c[0][1]->getSquareNorm() << '\n';
    std::cout << "Nabla_Psi_2c[0][2] = " << Nabla_Psi_2c[0][2]->getSquareNorm() << '\n';
    std::cout << "Nabla_Psi_2c[1][0] = " << Nabla_Psi_2c[1][0]->getSquareNorm() << '\n';
    std::cout << "Nabla_Psi_2c[1][1] = " << Nabla_Psi_2c[1][1]->getSquareNorm() << '\n';
    std::cout << "Nabla_Psi_2c[1][2] = " << Nabla_Psi_2c[1][2]->getSquareNorm() << '\n';
    // Nabla_Psi_2c[0][0] = Nabla_Psi_t[0]
    }

    

    std::vector<mrcpp::CompFunction<3> *> Nabla_Psi_t = Nabla_Psi_2c[0];
    std::vector<mrcpp::CompFunction<3> *> Nabla_Psi_b = Nabla_Psi_2c[1];


        

    // Compute Psi_top * K
    mrcpp::multiply<3>(Psi_t_K, Psi_2c[0], K_tree,building_precision);
    // Compute Psi_bottom * K
    mrcpp::multiply<3>(Psi_b_K, Psi_2c[1], K_tree,building_precision);

      

    if (debug){
        std:: cout << "Psi_t_K = " << Psi_t_K.getSquareNorm() << '\n';
        std:: cout << "Psi_b_K = " << Psi_b_K.getSquareNorm() << '\n';

        std::cout << '\n';
        std::cout << "Psi_t_K = " << Psi_t_K.getSquareNorm() << '\n';
        std::cout << "Psi_b_K = " << Psi_b_K.getSquareNorm() << '\n' << '\n';
        std::cout << "Nabla_Psi_t[0] = " << Nabla_Psi_t[0]->getSquareNorm() << '\n';
        std::cout << "Nabla_Psi_t[1] = " << Nabla_Psi_t[1]->getSquareNorm() << '\n';
        std::cout << "Nabla_Psi_t[2] = " << Nabla_Psi_t[2]->getSquareNorm() << '\n';
        std::cout << "Nabla_Psi_b[0] = " << Nabla_Psi_b[0]->getSquareNorm() << '\n';
        std::cout << "Nabla_Psi_b[1] = " << Nabla_Psi_b[1]->getSquareNorm() << '\n';
        std::cout << "Nabla_Psi_b[2] = " << Nabla_Psi_b[2]->getSquareNorm() << '\n';
        std::cout << '\n';
        std::cout << "Nabla_Psi_2c[0][0] = " << Nabla_Psi_2c[0][0]->getSquareNorm() << '\n';
        std::cout << "Nabla_Psi_2c[0][1] = " << Nabla_Psi_2c[0][1]->getSquareNorm() << '\n';
        std::cout << "Nabla_Psi_2c[0][2] = " << Nabla_Psi_2c[0][2]->getSquareNorm() << '\n';
        std::cout << "Nabla_Psi_2c[1][0] = " << Nabla_Psi_2c[1][0]->getSquareNorm() << '\n';
        std::cout << "Nabla_Psi_2c[1][1] = " << Nabla_Psi_2c[1][1]->getSquareNorm() << '\n';
        std::cout << "Nabla_Psi_2c[1][2] = " << Nabla_Psi_2c[1][2]->getSquareNorm() << '\n';
        std::cout << '\n';
        std::cout << '\n';
    }

    // Operatr ABGV
    mrcpp::ABGVOperator<3> D(MRA, 0.0, 0.0);

    // Gradient of | Psi_top * K >
    auto Nabla_Psi_t_K = mrcpp::gradient(D,Psi_t_K);
    // Gradient of Psi_bottom * K
    auto Nabla_Psi_b_K = mrcpp::gradient(D,Psi_b_K);

    // < \Nabla (\Psi * K) | \Nabla \Psi > 
    ComplexDouble Top_contribute = dot(*Nabla_Psi_t_K[0],*Nabla_Psi_t[0]) + dot(*Nabla_Psi_t_K[1],*Nabla_Psi_t[1]) + dot(*Nabla_Psi_t_K[2],*Nabla_Psi_t[2]);
    ComplexDouble Bottom_contribute = dot(*Nabla_Psi_b_K[0],*Nabla_Psi_b[0]) + dot(*Nabla_Psi_b_K[1],*Nabla_Psi_b[1]) + dot(*Nabla_Psi_b_K[2],*Nabla_Psi_b[2]);
    
    if (debug){
        std::cout << "Top_contribute = " << Top_contribute << '\n';
        std::cout << "Bottom_contribute = " << Bottom_contribute << '\n';
    }
    // REMEMBER THE MINUS SIGN!!! (it comes from the integration by parts)
    return -(Top_contribute + Bottom_contribute);
}



ComplexDouble compute_Term2_T_ZORA(MultiResolutionAnalysis<3> &MRA, std::vector<std::vector<mrcpp::CompFunction<3>*>> &Nabla_Psi_2c, mrcpp::CompFunction<3> &K_tree, std::vector<mrcpp::CompFunction<3> *> &Nabla_K_tree,  std::vector<mrcpp::CompFunction<3>> &Psi_2c){
    /// NEW_APPROACH: I first do all the \nabla(K) * \nabla(\Psi) and then multiply by the Psi from the left as a bra
    std::vector<mrcpp::CompFunction<3>> Nabla_K_Nabla_Psi_top(3, MRA);
    std::vector<mrcpp::CompFunction<3>> Nabla_K_Nabla_Psi_bottom(3, MRA);

    // \Naabla(K) * |\Nabla(\Psi)>
    /*
    for (int i=0;i<3;i++){
        mrcpp::multiply(Nabla_K_Nabla_Psi_top[i], *Nabla_K_tree[i], *Nabla_Psi_2c[0][i],building_precision, false, false, false);   // <---- Here we dereference the pointer
    }
    for (int i=0;i<3;i++){
        mrcpp::multiply(Nabla_K_Nabla_Psi_bottom[i], *Nabla_K_tree[i], *Nabla_Psi_2c[1][i],building_precision, false, false, false);   // <---- Here we dereference the pointer
    }

    */

    CompFunction<3> Nabla_K_Nabla_Psi_XTOP(MRA);
    CompFunction<3> Nabla_K_Nabla_Psi_YTOP(MRA);
    CompFunction<3> Nabla_K_Nabla_Psi_ZTOP(MRA);
    CompFunction<3> Nabla_K_Nabla_Psi_XBOTTOM(MRA);
    CompFunction<3> Nabla_K_Nabla_Psi_YBOTTOM(MRA);
    CompFunction<3> Nabla_K_Nabla_Psi_ZBOTTOM(MRA);

    // Let's do explicitely what the loops do:
    mrcpp::multiply(Nabla_K_Nabla_Psi_XTOP, *Nabla_K_tree[0], *Nabla_Psi_2c[0][0],building_precision, false, false, false);   
    mrcpp::multiply(Nabla_K_Nabla_Psi_YTOP, *Nabla_K_tree[1], *Nabla_Psi_2c[0][1],building_precision, false, false, false);
    mrcpp::multiply(Nabla_K_Nabla_Psi_ZTOP, *Nabla_K_tree[2], *Nabla_Psi_2c[0][2],building_precision, false, false, false);
    mrcpp::multiply(Nabla_K_Nabla_Psi_XBOTTOM, *Nabla_K_tree[0], *Nabla_Psi_2c[1][0],building_precision, false, false, false);
    mrcpp::multiply(Nabla_K_Nabla_Psi_YBOTTOM, *Nabla_K_tree[1], *Nabla_Psi_2c[1][1],building_precision, false, false, false);
    mrcpp::multiply(Nabla_K_Nabla_Psi_ZBOTTOM, *Nabla_K_tree[2], *Nabla_Psi_2c[1][2],building_precision, false, false, false);

    


    // Debug
    if (debug){
        std::cout << '\n';
        std::cout << "Nabla_K_Nabla_Psi_top[0] = " << Nabla_K_Nabla_Psi_top[0].getSquareNorm() << '\n';
        std::cout << "Nabla_K_Nabla_Psi_top[1] = " << Nabla_K_Nabla_Psi_top[1].getSquareNorm() << '\n';
        std::cout << "Nabla_K_Nabla_Psi_top[2] = " << Nabla_K_Nabla_Psi_top[2].getSquareNorm() << '\n';
        std::cout << "Nabla_K_Nabla_Psi_bottom[0] = " << Nabla_K_Nabla_Psi_bottom[0].getSquareNorm() << '\n';
        std::cout << "Nabla_K_Nabla_Psi_bottom[1] = " << Nabla_K_Nabla_Psi_bottom[1].getSquareNorm() << '\n';
        std::cout << "Nabla_K_Nabla_Psi_bottom[2] = " << Nabla_K_Nabla_Psi_bottom[2].getSquareNorm() << '\n';
        std::cout << '\n';
    }


    //ComplexDouble Top_contribution  = dot(Psi_2c[0],Nabla_K_Nabla_Psi_top[0]) + dot(Psi_2c[0],Nabla_K_Nabla_Psi_top[1]) + dot(Psi_2c[0],Nabla_K_Nabla_Psi_top[2]);
    //ComplexDouble Bottom_contribution  = dot(Psi_2c[1],Nabla_K_Nabla_Psi_bottom[0]) + dot(Psi_2c[1],Nabla_K_Nabla_Psi_bottom[1]) + dot(Psi_2c[1],Nabla_K_Nabla_Psi_bottom[2]);

    ComplexDouble Top_contribution  = dot(Psi_2c[0],Nabla_K_Nabla_Psi_XTOP) + dot(Psi_2c[0],Nabla_K_Nabla_Psi_YTOP) + dot(Psi_2c[0],Nabla_K_Nabla_Psi_ZTOP);
    ComplexDouble Bottom_contribution  = dot(Psi_2c[1],Nabla_K_Nabla_Psi_XBOTTOM) + dot(Psi_2c[1],Nabla_K_Nabla_Psi_YBOTTOM) + dot(Psi_2c[1],Nabla_K_Nabla_Psi_ZBOTTOM);


    return Top_contribution + Bottom_contribution;
}

void Cross_Product_Compfunct( std::vector<mrcpp::CompFunction<3>*> &Cross , std::vector<mrcpp::CompFunction<3>*> &A, std::vector<mrcpp::CompFunction<3>*> &B, MultiResolutionAnalysis<3> &MRA){
    Eigen::Matrix<int, 3, 2> Cross_coef;
    //          [A][B]
    Cross_coef << 1, 2,
                 2, 0,
                 0, 1;
                      
    CompFunction<3> Term_1(MRA);
    CompFunction<3> Term_2(MRA);

    if (debug){
        std::cout << "Cross_Product_Compfunct" << '\n';
        std::cout << "A[0] = " << A[0]->getSquareNorm() << '\n';
        std::cout << "A[1] = " << A[1]->getSquareNorm() << '\n';
        std::cout << "A[2] = " << A[2]->getSquareNorm() << '\n';
        std::cout << "B[0] = " << B[0]->getSquareNorm() << '\n';
        std::cout << "B[1] = " << B[1]->getSquareNorm() << '\n';
        std::cout << "B[2] = " << B[2]->getSquareNorm() << '\n';    
    }


    for (int i=0; i<3; i++){
        mrcpp::multiply(Term_1, *A[Cross_coef(i,0)], *B[Cross_coef(i,1)],building_precision, false, false, false); // A_i * B_j
        mrcpp::multiply(Term_2, *A[Cross_coef(i,1)], *B[Cross_coef(i,0)],building_precision, false, false, false); // A_j * B_i
        mrcpp::add<3>(*Cross[i], 1.0, Term_1 , -1.0, Term_2, building_precision, false); // \dv{F_A}{B} - \dv{F_B}{A}
    }

}




void compute_sigma_cdot_spinor(MultiResolutionAnalysis<3> &MRA  , std::vector<mrcpp::CompFunction<3>*> &Input_TOP, std::vector<mrcpp::CompFunction<3>*> &Input_BOTTOM, CompFunction<3> &SOC_Psi_t, CompFunction<3> &SOC_Psi_b){
    // In general, for each matrix vector(spinor) multiplication we need 4 scalar multiplication and 2 additions
    /*
    * SOC = sigma \cdot \binom{\vec{T}}{\vec{B}} 
    * This results in a 2 component spinor with the following operations:
    * [SOC TOP] = [B_x - i * B_y + T_z]
    * [SOC BOT] = [T_x + i * T_y - B_z]
    * 
    */

    std::vector<mrcpp::CompFunction<3>> SOC_yz(2); // This is the sum of the y and z components, temporary, to be used in the final sum
    //std::vector<mrcpp::CompFunction<3>> SOC_yz(2,MRA); see the difference? Well this doesn't work, and the one above yes
    // I'm not surprised that C++ programmers are so angry with the language, it is a mess
    // C++ is a language that is not meant to be used by humans, it is meant to be used by machines, and it is a mess
    // I don't know why I am using it, but I have to, because it is the only language that is used in this project
    // I don't like it, but I have to use it, so I will use it, and I will try to make it work
    // I will try to make it work, and I will try to make it work in a way that is understandable by humans
    // Even Copilot is confused by this, it doesn't know what to do with this code, it is a mess


    // A good recipe for a pasta alla carbonara is:
    // 1. Boil the pasta in salted water
    // 2. In a pan, cook the guanciale until crispy
    // 3. In a bowl, whisk the eggs with the cheese and pepper
    // 4. When the pasta is al dente, drain it and add it to the pan with the guanciale
    // 5. Remove the pan from the heat and add the egg mixture, stirring quickly to avoid scrambling the eggs
    // 6. Serve immediately with more cheese and pepper on top


    // Thus:
    // SOC_yz TOP = - i * B_y + T_z
    // SOC_yz BOT = + i * T_y - B_z

    mrcpp::add(SOC_yz[0], ComplexDouble(0.,-1.) , *Input_BOTTOM[1] , +1.0, *Input_TOP[2]    , building_precision, false);
    mrcpp::add(SOC_yz[1], ComplexDouble(0.,+1.) , *Input_TOP[1]    , -1.0, *Input_BOTTOM[2] , building_precision, false);

    // Finally i sum the x component to the yz component to the total SOC
    // SOC TOP = B_x + SOC_yz TOP
    // SOC BOT = T_x + SOC_yz BOT
    mrcpp::add(SOC_Psi_t, +1.0, *Input_BOTTOM[0], +1.0, SOC_yz[0],building_precision, false);
    mrcpp::add(SOC_Psi_b, +1.0, *Input_TOP[0], +1.0, SOC_yz[1],building_precision, false);

    if (debug){
        std::cout << "> NOW COMPUTING THE SIGMA CDOT SPINOR" << '\n';
        std::cout << '\n';
        std::cout << "Input_TOP[0] = " << Input_TOP[0]->getSquareNorm() << '\n';
        std::cout << "Input_TOP[1] = " << Input_TOP[1]->getSquareNorm() << '\n';
        std::cout << "Input_TOP[2] = " << Input_TOP[2]->getSquareNorm() << '\n';
        std::cout << "Input_BOTTOM[0] = " << Input_BOTTOM[0]->getSquareNorm() << '\n';
        std::cout << "Input_BOTTOM[1] = " << Input_BOTTOM[1]->getSquareNorm() << '\n';
        std::cout << "Input_BOTTOM[2] = " << Input_BOTTOM[2]->getSquareNorm() << '\n';

        std::cout << '\n';
        std::cout << "SOC_yz    TOP = " << SOC_yz[0].getSquareNorm() << '\n';
        std::cout << "SOC_yz BOTTOM = " << SOC_yz[1].getSquareNorm() << '\n';
        std::cout << '\n';
        std::cout << "SOC_Psi_t = " << SOC_Psi_t.getSquareNorm() << '\n';
        std::cout << "SOC_Psi_b = " << SOC_Psi_b.getSquareNorm() << '\n';
    }

}





ComplexDouble compute_Term3_T_ZORA(MultiResolutionAnalysis<3> &MRA, std::vector<std::vector<mrcpp::CompFunction<3>*>> &Nabla_Psi_2c, mrcpp::CompFunction<3> &K_tree, std::vector<mrcpp::CompFunction<3> *> &Nabla_K_tree,  std::vector<mrcpp::CompFunction<3>> &Psi_2c){
    // Compute the product K * (\Nabla \Psi) for top and bottom components
    std::vector<mrcpp::CompFunction<3>*> K_Nabla_Psi_top(3, new mrcpp::CompFunction<3>(MRA));
    std::vector<mrcpp::CompFunction<3>*> K_Nabla_Psi_bottom(3, new mrcpp::CompFunction<3>(MRA));

    std::vector<mrcpp::CompFunction<3> *> Nabla_Psi_top = Nabla_Psi_2c[0];
    std::vector<mrcpp::CompFunction<3> *> Nabla_Psi_bottom = Nabla_Psi_2c[1];

    mrcpp::ABGVOperator<3> D(MRA, 0.0, 0.0);

    // Compute the scalar product between the rotor and the sigma matrix vector
    CompFunction<3> SOC_Psi_t(MRA);
    CompFunction<3> SOC_Psi_b(MRA);

    std::vector<mrcpp::CompFunction<3>*> cross_top(3, new mrcpp::CompFunction<3>(MRA));
    std::vector<mrcpp::CompFunction<3>*> cross_bottom(3, new mrcpp::CompFunction<3>(MRA));
    Cross_Product_Compfunct(cross_top, Nabla_K_tree, Nabla_Psi_top, MRA);
    Cross_Product_Compfunct(cross_bottom, Nabla_K_tree, Nabla_Psi_bottom, MRA);

    compute_sigma_cdot_spinor(MRA,cross_top , cross_bottom , SOC_Psi_t, SOC_Psi_b);
    

    // Now i do the innetr produc with the bra
    
    ComplexDouble Top_contribution = dot(Psi_2c[0], SOC_Psi_t);
    ComplexDouble Bottom_contribution = dot(Psi_2c[1], SOC_Psi_b);

    if (debug){

        std::cout << '\n';
        std::cout << "TERM 3 COMPUTATION" << '\n';
        std::cout << "Nabla_Psi_top[0] = " << Nabla_Psi_top[0]->getSquareNorm() << '\n';
        std::cout << "Nabla_Psi_top[1] = " << Nabla_Psi_top[1]->getSquareNorm() << '\n';
        std::cout << "Nabla_Psi_top[2] = " << Nabla_Psi_top[2]->getSquareNorm() << '\n';   
        std::cout << '\n';
        std::cout << "Nabla_Psi_bottom[0] = " << Nabla_Psi_bottom[0]->getSquareNorm() << '\n';
        std::cout << "Nabla_Psi_bottom[1] = " << Nabla_Psi_bottom[1]->getSquareNorm() << '\n';
        std::cout << "Nabla_Psi_bottom[2] = " << Nabla_Psi_bottom[2]->getSquareNorm() << '\n';
        std::cout << '\n';
        std::cout << "K_tree = " << K_tree.getSquareNorm() << '\n';
        std::cout << '\n';
        std::cout << "Nabla_K_tree[0] = " << Nabla_K_tree[0]->getSquareNorm() << '\n';
        std::cout << "Nabla_K_tree[1] = " << Nabla_K_tree[1]->getSquareNorm() << '\n';
        std::cout << "Nabla_K_tree[2] = " << Nabla_K_tree[2]->getSquareNorm() << '\n';
        std::cout << '\n';
        std::cout << "Cross top[0] = " << cross_top[0]->getSquareNorm() << '\t' << "Real: " << cross_top[0]->isreal() << "; Complex: " << cross_top[0]->iscomplex() << '\n';
        std::cout << "Cross top[1] = " << cross_top[1]->getSquareNorm() << '\t' << "Real: " << cross_top[1]->isreal() << "; Complex: " << cross_top[1]->iscomplex() << '\n';
        std::cout << "Cross top[2] = " << cross_top[2]->getSquareNorm() << '\t' << "Real: " << cross_top[2]->isreal() << "; Complex: " << cross_top[2]->iscomplex() << '\n';
        std::cout << "Cross bottom[0] = " << cross_bottom[0]->getSquareNorm() << '\t' << "Real: " << cross_bottom[0]->isreal() << "; Complex: " << cross_bottom[0]->iscomplex() << '\n';
        std::cout << "Cross bottom[1] = " << cross_bottom[1]->getSquareNorm() << '\t' << "Real: " << cross_bottom[1]->isreal() << "; Complex: " << cross_bottom[1]->iscomplex() << '\n';
        std::cout << "Cross bottom[2] = " << cross_bottom[2]->getSquareNorm() << '\t' << "Real: " << cross_bottom[2]->isreal() << "; Complex: " << cross_bottom[2]->iscomplex() << '\n';
        std::cout << '\n';
        std::cout << "SOC_Psi_t = " << SOC_Psi_t.getSquareNorm() << '\n';
        std::cout << "SOC_Psi_b = " << SOC_Psi_b.getSquareNorm() << '\n';
        std::cout << '\n';
        std::cout << "Top_contribution = " << Top_contribution << '\n';
        std::cout << "Bottom_contribution = " << Bottom_contribution << '\n';


        std::cout << '\n';

    }

    for (int i=0; i<3; i++){
        K_Nabla_Psi_top[i]->free(); // Clear the previous values
        //delete K_Nabla_Psi_top[i];  // Delete the object
        
        K_Nabla_Psi_bottom[i]->free(); // Clear the previous values
        //delete K_Nabla_Psi_bottom[i];  // Delete the object
        
        cross_top[i]->free(); // Clear the previous values
        //delete cross_top[i];  // Delete the object

        cross_bottom[i]->free(); // Clear the previous values
        //delete cross_bottom[i];  // Delete the object
    }
    K_Nabla_Psi_bottom.clear(); // Clear the outer container
    K_Nabla_Psi_top.clear(); // Clear the outer container
    cross_top.clear(); // Clear the outer container
    cross_bottom.clear(); // Clear the outer container
    



    return ComplexDouble(0,1)*(Top_contribution + Bottom_contribution);
}







void compute_energy_1e_ZORA(MultiResolutionAnalysis<3> &MRA, std::vector<std::vector<mrcpp::CompFunction<3>*>> &Nabla_Psi_2c, mrcpp::CompFunction<3> &Kappa_tree, std::vector<mrcpp::CompFunction<3> *> &Nabla_Kappa_tree,  Spinor_CompFunction &Psi_2c, CompFunction<3> &V, std::vector<mrcpp::CompFunction<3> *> J_Psi, std::vector<mrcpp::CompFunction<3> *> &K_Psi, double &Orbital_energy){
    if (verbose){std::cout << '\t' << "Computing Kinetic Term 1..." << '\n';}
    ComplexDouble Term1 = compute_Term1_T_ZORA(MRA, Nabla_Psi_2c, Kappa_tree, Psi_2c);
    if (verbose){std::cout << '\t' << "Computing Kinetic Term 2..." << '\n';}
    ComplexDouble Term2 = compute_Term2_T_ZORA(MRA, Nabla_Psi_2c, Kappa_tree, Nabla_Kappa_tree, Psi_2c);
    if (verbose){std::cout << '\t' << "Computing Kinetic Term 2..." << '\n';}
    ComplexDouble Term3 = compute_Term3_T_ZORA(MRA, Nabla_Psi_2c, Kappa_tree, Nabla_Kappa_tree, Psi_2c);

    // Now we do <Psi | V | Psi>
    mrcpp::CompFunction<3> psi_top__V(MRA);
    mrcpp::CompFunction<3> psi_bottom__V(MRA);

    // V*|Psi>
    mrcpp::multiply(psi_top__V, V, Psi_2c[0],building_precision, false, false, false);
    mrcpp::multiply(psi_bottom__V, V, Psi_2c[1],building_precision, false, false, false);

    if (verbose){std::cout << '\t' << "Computing Potential Term..." << '\n' <<'\n';}

    ComplexDouble monoel_pot_energy = dot(Psi_2c[0],psi_top__V) + dot(Psi_2c[1],psi_bottom__V);

    // BIELECTRONIC PART
    ComplexDouble Coulomb_repulsion = dot(Psi_2c[0],*J_Psi[0]) + dot(Psi_2c[1],*J_Psi[1]);
    ComplexDouble Exchange_attraction = dot(Psi_2c[0],*K_Psi[0]) + dot(Psi_2c[1],*K_Psi[1]);
    ComplexDouble biel_energy = 0.5 * (Coulomb_repulsion - Exchange_attraction);


    double kin_prefactor = -1.0/(2.0*m);
    ComplexDouble kinetic_energy = kin_prefactor * (Term1 + Term2 + Term3) ;

    double real_Z = double(Z);
    ComplexDouble total_energy =  kinetic_energy + monoel_pot_energy + biel_energy;

    verbose = true;
    if (verbose){
    std::cout << "--------------------------------------------------" << '\n';
    std::cout << "Kinetic energy = "  << kinetic_energy << '\n';
    std::cout << "  TERM 1 = "  << kin_prefactor * Term1 << '\n';
    std::cout << "  TERM 2 = "  << kin_prefactor * Term2 << '\n';
    std::cout << "  TERM 3 = "  << kin_prefactor * Term3 << '\n';
    
    std::cout << "Potential energy c-e  = "  << monoel_pot_energy << '\n';
    std::cout << "Potential energy e-e  = "  << biel_energy << '\n'<< '\n';

    std::cout << "Total  = "  << total_energy << '\n';
    std::cout << "--------------------------------------------------" << '\n';
    }

    verbose = false;    


    std::cout << '\t' << "-------------------------------" << '\n';
        std::cout << '\t' << "| Energy = " << total_energy.real() << "|" << '\n';
        std::cout << '\t' << "-------------------------------" << '\n' << '\n';


    // IF the enery is above zero it stops the program
    
    //if (total_energy.real() > 0){
    //   std::cout << "Energy is above zero, stopping the program" << '\n';
    //    exit(1);
    //}

    
    Orbital_energy = total_energy.real();

}


void compute_energy_atom_ZORA(MultiResolutionAnalysis<3> &MRA, Spinor_gradients_pointer_list &Nabla_Psi_List, mrcpp::CompFunction<3> &Kappa_tree, std::vector<mrcpp::CompFunction<3> *> &Nabla_Kappa_tree,  Orbital_Pointer_List &Psi_List, CompFunction<3> &V, std::vector<std::vector<mrcpp::CompFunction<3> *>> J_List, std::vector<std::vector<mrcpp::CompFunction<3> *>> &K_List, double &Atomic_energy){
    std::vector<double> Orbital_energy_list;
    double tmp;
    for (int i = 0; i<4 ; i++){
        compute_energy_1e_ZORA(MRA, Nabla_Psi_List[i], Kappa_tree, Nabla_Kappa_tree, *Psi_List[i], V, J_List[i], K_List[i], tmp);
        Orbital_energy_list.push_back(tmp);
        std::cout << '\t' << "Orbital energy " << i << " = " << tmp << '\n';
    }
    Atomic_energy = std::accumulate(Orbital_energy_list.begin(), Orbital_energy_list.end(), 0.0);
    std::cout << '\n';
    std::cout << '\t' << "Atomic energy = " << Atomic_energy << '\n';
    std::cout << '\n';
}


void compute_term_A_Propagator(MultiResolutionAnalysis<3> &MRA, double E_n, std::vector<mrcpp::CompFunction<3>> &Psi_in, CompFunction<3> &V,CompFunction<3> &term_A_top,CompFunction<3> &term_A_bottom){
    double const_factor = -E_n / (c*c);
    if (debug){
    std::cout << "--------------------------------------------------" << '\n';
    std::cout << " Input in A:" << '\n';
    std::cout << "  Psi_in[0] = " << Psi_in[0].getSquareNorm() << '\n';
    std::cout << "  Psi_in[1] = " << Psi_in[1].getSquareNorm() << '\n';
    std::cout << "  V = " << V.getSquareNorm() << '\n';
}

    mrcpp::multiply(building_precision, term_A_top, 1.0, V, Psi_in[0]);
    mrcpp::multiply(building_precision, term_A_bottom, 1.0, V, Psi_in[1]);
    term_A_top.rescale(const_factor);
    term_A_bottom.rescale(const_factor);

}


void compute_term_B_Propagator(MultiResolutionAnalysis<3> &MRA, std::vector<mrcpp::CompFunction<3>> &Psi_in, std::vector<std::vector<mrcpp::CompFunction<3> *>> &Nabla_Psi_2c, std::vector<mrcpp::CompFunction<3> *> &Nabla_K_tree , CompFunction<3> &K_tree, CompFunction<3> &K_inverse, CompFunction<3> &term_B_top,CompFunction<3> &term_B_bottom){
    std::vector<mrcpp::CompFunction<3>> Nabla_K_Nabla_Psi_top(3, MRA);
    std::vector<mrcpp::CompFunction<3>> Nabla_K_Nabla_Psi_bottom(3, MRA);

    // EXPLICIT RUTE
    mrcpp::CompFunction<3> DOT_TOP_x(MRA); // ----> X component for the dot product
    mrcpp::CompFunction<3> DOT_TOP_y(MRA); // ----> Y component for the dot product
    mrcpp::CompFunction<3> DOT_TOP_z(MRA); // ----> Z component for the dot product
    mrcpp::CompFunction<3> DOT_BOTTOM_x(MRA); // ----> X component for the dot product
    mrcpp::CompFunction<3> DOT_BOTTOM_y(MRA); // ----> Y component for the dot product
    mrcpp::CompFunction<3> DOT_BOTTOM_z(MRA); // ----> Z component for the dot product


    mrcpp::multiply(DOT_TOP_x, *Nabla_K_tree[0], *Nabla_Psi_2c[0][0],building_precision, false, false, false);
    mrcpp::multiply(DOT_TOP_y, *Nabla_K_tree[1], *Nabla_Psi_2c[0][1],building_precision, false, false, false);
    mrcpp::multiply(DOT_TOP_z, *Nabla_K_tree[2], *Nabla_Psi_2c[0][2],building_precision, false, false, false);
    mrcpp::multiply(DOT_BOTTOM_x, *Nabla_K_tree[0], *Nabla_Psi_2c[1][0],building_precision, false, false, false);
    mrcpp::multiply(DOT_BOTTOM_y, *Nabla_K_tree[1], *Nabla_Psi_2c[1][1],building_precision, false, false, false);
    mrcpp::multiply(DOT_BOTTOM_z, *Nabla_K_tree[2], *Nabla_Psi_2c[1][2],building_precision, false, false, false);

    // Now i sum the y and z components
    mrcpp::CompFunction<3> DOT_TOP_yz(MRA);
    mrcpp::CompFunction<3> DOT_BOTTOM_yz(MRA);
    mrcpp::CompFunction<3> DOT_TOP(MRA);
    mrcpp::CompFunction<3> DOT_BOTTOM(MRA);

    mrcpp::add(DOT_TOP_yz, 1.0, DOT_TOP_y, 1.0, DOT_TOP_z,building_precision, false);
    mrcpp::add(DOT_BOTTOM_yz, 1.0, DOT_BOTTOM_y, 1.0, DOT_BOTTOM_z,building_precision, false);

    // Finally i sum the x component to the yz component to the total SOC
    mrcpp::add(DOT_TOP, 1.0, DOT_TOP_x, 1.0, DOT_TOP_yz,building_precision, false);
    mrcpp::add(DOT_BOTTOM, 1.0, DOT_BOTTOM_x, 1.0, DOT_BOTTOM_yz,building_precision, false);

    // Multiply by K_inverse
    mrcpp::multiply(term_B_top, K_inverse, DOT_TOP,building_precision, false, false, false);
    mrcpp::multiply(term_B_bottom, K_inverse, DOT_BOTTOM,building_precision, false, false, false);
}


void compute_term_C_Propagator(MultiResolutionAnalysis<3> &MRA, std::vector<std::vector<mrcpp::CompFunction<3> *>> &Nabla_Psi_2c , std::vector<mrcpp::CompFunction<3> *> &Nabla_K_tree, CompFunction<3> &K_inverse, CompFunction<3> &term_C_top,CompFunction<3> &term_C_bottom){
    // Where I'll put the rotor of K * (\Nabla \Psi)

    // Now we compute the scalar product between the rotor and the sigma matrix vector
    CompFunction<3> SOC_Psi_t(MRA);
    CompFunction<3> SOC_Psi_b(MRA);
    
    std::vector<mrcpp::CompFunction<3>*> cross_top(3, new mrcpp::CompFunction<3>(MRA));
    std::vector<mrcpp::CompFunction<3>*> cross_bottom(3, new mrcpp::CompFunction<3>(MRA));
    
    std::vector<mrcpp::CompFunction<3> *> Nabla_Psi_top(3, new mrcpp::CompFunction<3>(MRA));
    std::vector<mrcpp::CompFunction<3> *> Nabla_Psi_bottom(3, new mrcpp::CompFunction<3>(MRA));

    for (int i=0;i<3;i++){
        Nabla_Psi_top[i] = Nabla_Psi_2c[0][i];
        Nabla_Psi_bottom[i] = Nabla_Psi_2c[1][i];
    }
 
    if (verbose){std::cout << '\t' << "Computing the cross product for the SOC (TOP)" << '\n';}
    Cross_Product_Compfunct(cross_top, Nabla_K_tree, Nabla_Psi_top, MRA);
    if (verbose){std::cout << '\t' << "Computing the cross product for the SOC (BOT)" << '\n' <<'\n';}
    Cross_Product_Compfunct(cross_bottom, Nabla_K_tree, Nabla_Psi_bottom, MRA);


    if (verbose){std::cout << '\t' << "Computing the dot product with the Pauli Matrices" << '\n' ;}
    compute_sigma_cdot_spinor(MRA, cross_top, cross_bottom, SOC_Psi_t, SOC_Psi_b);
    
    // Now the actual top and bottom components for C are the i * K_inverse * SOC_Psi
    // Times K INVERSE
    mrcpp::multiply(term_C_top, K_inverse, SOC_Psi_t,building_precision, false, false, false);
    mrcpp::multiply(term_C_bottom, K_inverse, SOC_Psi_b,building_precision, false, false, false);
    
    // Times i 
    if (debug){
    std::cout << "NORM PRE RESCALE:" << '\n';
    std::cout << "term_C_top = " << term_C_top.getSquareNorm() << '\n';
    std::cout << "term_C_bottom = " << term_C_bottom.getSquareNorm() << '\n';
    }

    for (int i=0; i<3; i++){
        Nabla_Psi_top[i]->free(); // Clear the previous values
        //delete K_Nabla_Psi_top[i];  // Delete the object
        
        Nabla_Psi_bottom[i]->free(); // Clear the previous values
        //delete K_Nabla_Psi_bottom[i];  // Delete the object
        
        cross_top[i]->free(); // Clear the previous values
        //delete cross_top[i];  // Delete the object

        cross_bottom[i]->free(); // Clear the previous values
        //delete cross_bottom[i];  // Delete the object
    }
    Nabla_Psi_bottom.clear(); // Clear the outer container
    Nabla_Psi_top.clear(); // Clear the outer container
    cross_top.clear(); // Clear the outer container
    cross_bottom.clear(); // Clear the outer container
   



    term_C_top.rescale(ComplexDouble(0,1));
    term_C_bottom.rescale(ComplexDouble(0,1));
    
}


void compute_term_D_Propagator(MultiResolutionAnalysis<3> &MRA,std::vector<mrcpp::CompFunction<3>> &Psi_in ,CompFunction<3> &V , CompFunction<3> &K_inverse, CompFunction<3> &term_D_top,CompFunction<3> &term_D_bottom){
    // Term_D = (-2m/K) * V * Psi
    double factor = -2.0 * m;
    CompFunction<3> V_Psi_t(MRA);
    CompFunction<3> V_Psi_b(MRA);

    
    mrcpp::multiply(building_precision, V_Psi_t, 1.0, V, Psi_in[0]);
    mrcpp::multiply(building_precision, V_Psi_b, 1.0, V, Psi_in[1]);
    
    if (debug){
    std::cout << "NORM INPUT:" << '\n';
    std::cout << "V_Psi_top = " << V_Psi_t.getSquareNorm() << '\n';
    std::cout << "V_Psi_bottom = " << V_Psi_b.getSquareNorm() << '\n';
    }
    
    mrcpp::multiply(building_precision, term_D_top ,    1.0, K_inverse,  V_Psi_t);
    mrcpp::multiply(building_precision, term_D_bottom , 1.0, K_inverse, V_Psi_b);



    // I don't know Walt... kinda sus
    term_D_top.rescale(factor);
    term_D_bottom.rescale(factor);
}


void compute_term_E_Propagator(MultiResolutionAnalysis<3> &MRA,std::vector<std::vector<std::complex<double>>> &F_matrix, int Orbital_Index, Orbital_Pointer_List Psi_List, CompFunction<3> &K_inverse, CompFunction<3> &term_E_top,CompFunction<3> &term_E_bottom){

    std::vector<mrcpp::CompFunction<3>> Psi_tops;
    std::vector<mrcpp::CompFunction<3>> Psi_bottoms;
    std::vector<ComplexDouble> coeffs;

    for (int i=0; i<3; i++){
        if (i== Orbital_Index){
            continue; // Skip the orbital we are computing the energy for
        }
        Psi_tops.push_back((*Psi_List[i])[0]);
        Psi_bottoms.push_back((*Psi_List[i])[1]);
        coeffs.push_back(2*m*F_matrix[Orbital_Index][i]); // We want to sum all the other orbitals
    }

    CompFunction<3> Psi_top_sum(MRA);
    CompFunction<3> Psi_bottom_sum(MRA);
    
    linear_combination(Psi_top_sum, coeffs, Psi_tops, building_precision, false);
    linear_combination(Psi_bottom_sum, coeffs, Psi_bottoms, building_precision, false);


    mrcpp::multiply(term_E_top, K_inverse, Psi_top_sum, building_precision, false, false, false);
    mrcpp::multiply(term_E_bottom, K_inverse, Psi_bottom_sum, building_precision, false, false, false);

    
}



void apply_Helmholtz_ZORA(MultiResolutionAnalysis<3> &MRA,
                          std::vector<std::vector<std::complex<double>>> &F_matrix, int Orbital_Index,
                          Orbital_Pointer_List &Psi_list, std::vector<mrcpp::CompFunction<3>> &Psi_in, 
                          std::vector<std::vector<mrcpp::CompFunction<3> *>> &Nabla_Psi_2c, 
                          CompFunction<3> &V, 
                          Spinor_CompFunction J_Psi_top,
                          Spinor_CompFunction J_Psi_bottom,
                          Spinor_CompFunction K_Psi_top,
                          Spinor_CompFunction K_Psi_bottom,
                          std::vector<mrcpp::CompFunction<3> *> &Nabla_K_tree,  
                          CompFunction<3> &K_tree, 
                          CompFunction<3> &K_inverse, 
                          double &E_n, 
                          std::vector<mrcpp::CompFunction<3>> &Psi_out){
    // Build the Helmholtz operator
    
    double mu = std::sqrt(-2*E_n);
    if (verbose || debug){
    std::cout << "--------------------------------------------------" << '\n';
    std::cout << "Mu = " << mu << '\n';
    std::cout << "--------------------------------------------------" << '\n';
    }
    mrcpp::HelmholtzOperator Helm(MRA, mu, building_precision);

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

    // Term a) -E^n/c^2 * V * Psi^n
    CompFunction<3> term_A_top(MRA);
    CompFunction<3> term_A_bottom(MRA);

    // Term b) + K_invinverted * (\Nabla K \cdot \Nabla \Psi^n)
    CompFunction<3> term_B_top(MRA);
    CompFunction<3> term_B_bottom(MRA);

    // Term c) + K_invinverted * (i * \sigma \cdot (\Nabla \cross K (\Nabla (\Psi^n))))
    CompFunction<3> term_C_top(MRA);
    CompFunction<3> term_C_bottom(MRA);

    // Term d) - K_invinverted * 2*m*V
    CompFunction<3> term_D_top(MRA);
    CompFunction<3> term_D_bottom(MRA);
    CompFunction<3> term_E_top(MRA);
    CompFunction<3> term_E_bottom(MRA);

    // Compute the term A
    if (verbose){std::cout << "Computing term A..." << '\n';}
    compute_term_A_Propagator(MRA, E_n, Psi_in, V, term_A_top, term_A_bottom);

    // Compute the term B
    if (verbose){std::cout << "Computing term B..." << '\n';}
    compute_term_B_Propagator(MRA, Psi_in, Nabla_Psi_2c, Nabla_K_tree, K_tree, K_inverse, term_B_top, term_B_bottom);
    
    // Compute the term C
    if (verbose){std::cout << "Computing term C..." << '\n';}
    compute_term_C_Propagator(MRA, Nabla_Psi_2c, Nabla_K_tree, K_inverse, term_C_top, term_C_bottom);

    // Compute the term D
    if (verbose){std::cout << "Computing term D..." << '\n';}
    compute_term_D_Propagator(MRA, Psi_in, V, K_inverse, term_D_top, term_D_bottom);

    // Compute the term E
    if (verbose){std::cout << "Computing term E..." << '\n';}
    compute_term_E_Propagator(MRA, F_matrix, Orbital_Index, Psi_list, K_inverse, term_E_top, term_E_bottom);

    // Now we sum all the terms, we have a total of 4 terms to sum. We'll take them by pairs and sum them
    CompFunction<3> first_pair_top(MRA);
    CompFunction<3> first_pair_bottom(MRA);
    CompFunction<3> second_pair_top(MRA);
    CompFunction<3> second_pair_bottom(MRA);

    // DEBUG
    if (debug){
    std::cout << '\n' << "Norm of the components" << '\n';
    std::cout << "term_A_top = " << term_A_top.getSquareNorm() << '\t' << "term_A_bottom = " << term_A_bottom.getSquareNorm() << '\n';
    std::cout << "term_B_top = " << term_B_top.getSquareNorm() << '\t' << "term_B_bottom = " << term_B_bottom.getSquareNorm() << '\n';
    std::cout << "term_C_top = " << term_C_top.getSquareNorm() << '\t' << "term_C_bottom = " << term_C_bottom.getSquareNorm() << '\n';
    std::cout << "term_D_top = " << term_D_top.getSquareNorm() << '\t' << "term_D_bottom = " << term_D_bottom.getSquareNorm() << '\n';
    std::cout << "term_E_top = " << term_E_top.getSquareNorm() << '\t' << "term_E_bottom = " << term_E_bottom.getSquareNorm() << '\n';
    }


    // Adding all together
    if (verbose){std::cout << '\n'<< "Adding up all the contributions..." << '\n';}
    mrcpp::add(first_pair_top, 1.0, term_A_top, 1.0, term_B_top, building_precision, false);
    mrcpp::add(first_pair_bottom, 1.0, term_A_bottom, 1.0, term_B_bottom, building_precision, false);
    mrcpp::add(second_pair_top, 1.0, term_C_top, 1.0, term_D_top, building_precision, false);
    mrcpp::add(second_pair_bottom, 1.0, term_C_bottom, 1.0, term_D_bottom, building_precision, false);

    // SECOND DEBUG
    if (debug){
    std::cout << '\n' << "Norm of the components after the addition" << '\n';
    std::cout << "first_pair_top = " << first_pair_top.getSquareNorm() << '\t' << "first_pair_bottom = " << first_pair_bottom.getSquareNorm() << '\n';
    std::cout << "second_pair_top = " << second_pair_top.getSquareNorm() << '\t' << "second_pair_bottom = " << second_pair_bottom.getSquareNorm() << '\n';
    }

    // Now we sum the two pairs
    std::vector<mrcpp::CompFunction<3>> Psi_to_be_convoluted(2,MRA);

    // Finalizing the sum of all terms, i put them in 2 temporary CompFunction
    mrcpp::CompFunction<3> tmp_top(MRA);
    mrcpp::CompFunction<3> tmp_bot(MRA);
    mrcpp::add(tmp_top, 1.0, first_pair_top, 1.0, second_pair_top, building_precision, false);
    mrcpp::add(tmp_bot, 1.0, first_pair_bottom, 1.0, second_pair_bottom, building_precision, false);


    mrcpp::CompFunction<3> tmp_top_2(MRA);
    mrcpp::CompFunction<3> tmp_bot_2(MRA);
    mrcpp::add(tmp_top_2, 1.0, tmp_top, 1.0, term_E_top, building_precision, false);
    mrcpp::add(tmp_bot_2, 1.0, tmp_bot, 1.0, term_E_bottom, building_precision, false);

    // printing the norm of top and bottom compoennts for  tmp an term_E
    
    if (debug){
    std::cout << '\n' << "Norm of the components after the final addition" << '\n';
    std::cout << "tmp_top = " << tmp_top.getSquareNorm() << '\t' << "tmp_bot = " << tmp_bot.getSquareNorm() << '\n';
    std::cout << "term_E_top = " << term_E_top.getSquareNorm() << '\t' << "term_E_bottom = " << term_E_bottom.getSquareNorm() << '\n';
    std::cout << "tmp_top_2 = " << tmp_top_2.getSquareNorm() << '\t' << "tmp_bot_2 = " << tmp_bot_2.getSquareNorm() << '\n';
    }


    if (verbose){std::cout << '\n' << "Finalizing the bielectronic contributions..." << '\n';}
    // J-K|Psi> part
    CompFunction<3> bielectronic_term_top(MRA);
    CompFunction<3> bielectronic_term_bottom(MRA);

    // (J-K)|Psi[i]>
    mrcpp::add(bielectronic_term_top, 1.0, J_Psi_top[Orbital_Index], -1.0, K_Psi_top[Orbital_Index], building_precision, false);
    mrcpp::add(bielectronic_term_bottom, 1.0, J_Psi_bottom[Orbital_Index], -1.0, K_Psi_bottom[Orbital_Index], building_precision, false);

    // -2m * (J-K)|Psi[i]>
    bielectronic_term_top.rescale(ComplexDouble(-2*m,0));
    bielectronic_term_bottom.rescale(ComplexDouble(-2*m,0));


    // Times kappa^{-1}
    CompFunction<3> bielectronic_term_K_inverse_top(MRA);
    CompFunction<3> bielectronic_term_K_inverse_bottom(MRA);
    
    // -2m * (J-K)|Psi[i]> / Kappa
    mrcpp::multiply(bielectronic_term_K_inverse_top, K_inverse, bielectronic_term_top, building_precision, false, false, false);
    mrcpp::multiply(bielectronic_term_K_inverse_bottom, K_inverse, bielectronic_term_bottom, building_precision, false, false, false);


    //Print the norms for the bielectronic terms
    
    if (debug){
    std::cout << '\n' << "Norm of the bielectronic terms" << '\n';
    std::cout << "bielectronic_term_top = " << bielectronic_term_top.getSquareNorm() << '\n';
    std::cout << "bielectronic_term_bottom = " << bielectronic_term_bottom.getSquareNorm() << '\n';
    std::cout << "bielectronic_term_K_inverse_top = " << bielectronic_term_K_inverse_top.getSquareNorm() << '\n';
    std::cout << "bielectronic_term_K_inverse_bottom = " << bielectronic_term_K_inverse_bottom.getSquareNorm() << '\n';
    }



    // Add the bielectronic term to the rest of the contributions in Psi_to_be_convoluted
    if (verbose){std::cout << '\n' << "Adding the bielectronic contributions..." << '\n';}
    CompFunction<3> Psi_to_be_convoluted_top(MRA);
    CompFunction<3> Psi_to_be_convoluted_bottom(MRA);
    mrcpp::add(Psi_to_be_convoluted_top, 1.0, tmp_top_2, 1.0, bielectronic_term_K_inverse_top, building_precision, false);
    mrcpp::add(Psi_to_be_convoluted_bottom, 1.0, tmp_bot_2, 1.0, bielectronic_term_K_inverse_bottom, building_precision, false);


    // LAST DEBUG
    if (debug){
    std::cout << '\n' << "Norm of the components after the final addition" << '\n';
    std::cout << "Psi_to_be_convoluted[0] = " << Psi_to_be_convoluted[0].getSquareNorm() << '\n' << "Psi_to_be_convoluted[1] = " << Psi_to_be_convoluted[1].getSquareNorm() << '\n';
    }


    // Now we apply the Helmholtz operator
    if (verbose){std::cout << '\n' <<  "!APPLYING HELMOLTZ OPERATOR!" << '\n' << '\n';}
    CompFunction<3> Psi_out_top(MRA);
    CompFunction<3> Psi_out_bottom(MRA);

    mrcpp::apply(building_precision, Psi_out_top, Helm, Psi_to_be_convoluted_top);
    mrcpp::apply(building_precision, Psi_out_bottom, Helm, Psi_to_be_convoluted_bottom);
    std::cout << "Applyed!" << '\n';

    std::cout << '\n';
    // print the norm of the output components
    if (debug){
    std::cout << "Norm of the output components:" << '\n';
    std::cout << "Psi_out[0] = " << Psi_out_top.getSquareNorm() << '\n';
    std::cout << "Psi_out[1] = " << Psi_out_bottom.getSquareNorm() << '\n';
    }

    Psi_out[0] = Psi_out_top;
    Psi_out[1] = Psi_out_bottom;

}


void apply_Helmholtz_ZORA_all_electrons(MultiResolutionAnalysis<3> &MRA, 
    std::vector<std::vector<std::complex<double>>> &F_matrix, 
    Orbital_Pointer_List &Psi_list, 
    Spinor_gradients_pointer_list Nabla_psi_list, 
    CompFunction<3> &V, 
    Spinor_CompFunction J_Psi_top,
    Spinor_CompFunction J_Psi_bottom,
    Spinor_CompFunction K_Psi_top,
    Spinor_CompFunction K_Psi_bottom,
    std::vector<mrcpp::CompFunction<3> *> &Nabla_K_tree,  
    CompFunction<3> &K_tree, 
    CompFunction<3> &K_inverse,  
    std::vector<double> &norm_list){

    // This function is the same as the previous one, but it takes a Orbital_Pointer_List as input and output
    // It is used to apply the Helmholtz operator to the orbitals of the atom
    //Orbital_Pointer_List Psi_out_vector(4);
    Spinor_CompFunction Psi_out_i(2,MRA);
    std::vector<CompFunction<3>> Psi_out_list_top(4);
    std::vector<CompFunction<3>> Psi_out_list_bottom(4);

    Spinor_CompFunction Psi_in(2, MRA);
    CompFunction<3> difference_top(MRA);
    CompFunction<3> difference_bottom(MRA);

    double energy_fock = 0.0;

    if (verbose){std::cout << "Applying Helmholtz operator to all electrons..." << '\n';}
    
     for (int i = 0; i < 4; i++){
        std::cout << "-------------------------------------------" << '\n';
        std::cout << "Applying Helmholtz operator to orbital " << i << "..." << '\n';
        Psi_in[0] = (*Psi_list[i])[0];
        Psi_in[1] = (*Psi_list[i])[1];
        energy_fock = F_matrix[i][i].real();
        apply_Helmholtz_ZORA(MRA, F_matrix, i, Psi_list, Psi_in, Nabla_psi_list[i],V, J_Psi_top, J_Psi_bottom, K_Psi_top, K_Psi_bottom, Nabla_K_tree, K_tree, K_inverse, energy_fock, Psi_out_i);
        if (debug){
        std::cout << "Norm of the output components:" << '\n';
        std::cout << "Psi_out_vector[" << i << "][0] = " << Psi_out_i[0].getSquareNorm() << '\n';
        std::cout << "Psi_out_vector[" << i << "][1] = " << Psi_out_i[1].getSquareNorm() << '\n';
        std::cout << "-------------------------------------------" << '\n';
        // also print the adresses of the output components
        std::cout << "Psi_out_vector[" << i << "][0] address = " << &Psi_out_i[0] << '\n';
        std::cout << "Psi_out_vector[" << i << "][1] address = " << &Psi_out_i[1] << '\n';
        std::cout << "-------------------------------------------" << '\n';
        }
        


        // Now we compute the difference between the output and the input
        mrcpp::add(difference_top, 1.0, Psi_out_i[0], -1.0, (*Psi_list[i])[0], building_precision, false);
        mrcpp::add(difference_bottom, 1.0, Psi_out_i[1], -1.0, (*Psi_list[i])[1], building_precision, false);
       

        norm_list[i] = std::sqrt(difference_top.getSquareNorm() + difference_bottom.getSquareNorm());
        // Store in a temporary array the now modified convoluted spinor functions
        deep_copy(Psi_out_list_top[i] , Psi_out_i[0]);
        deep_copy(Psi_out_list_bottom[i], Psi_out_i[1]);
    }

    // Once we have applied the Helmholtz operator to all the orbitals, we can copy the convoluted results in the actual Psi_list
    for (int i = 0; i < 4; i++){
        deep_copy((*Psi_list[i])[0], Psi_out_list_top[i]);
        deep_copy((*Psi_list[i])[1], Psi_out_list_bottom[i]);
    }

}

void Renormalize_Spinor(std::vector<mrcpp::CompFunction<3>> &Psi_2c){
    double norm = std::sqrt(Psi_2c[0].getSquareNorm() + Psi_2c[1].getSquareNorm());
    Psi_2c[0].rescale(1.0/norm);
    Psi_2c[1].rescale(1.0/norm);

    norm = std::sqrt(Psi_2c[0].getSquareNorm() + Psi_2c[1].getSquareNorm());
    if (norm < 0.999 || norm > 1.001){
        std::cerr << "ERROR : Norm is not 1.0, it is " << norm << '\n';
    }
}





















