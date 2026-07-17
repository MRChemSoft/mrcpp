

double compute_1e_Kinetic_energy(MultiResolutionAnalysis<3> &MRA, FunctionTree<3> &Phi){
    mrcpp::ABGVOperator<3> D(MRA, 0.0, 0.0);
    auto Nabla_Phi = mrcpp::gradient(D, Phi);



    double x_component = mrcpp::dot(get_func(Nabla_Phi,0), get_func(Nabla_Phi,0));
    double y_component = mrcpp::dot(get_func(Nabla_Phi,1), get_func(Nabla_Phi,1));
    double z_component = mrcpp::dot(get_func(Nabla_Phi,2), get_func(Nabla_Phi,2));

    double kinetic_energy = 0.5 * (x_component + y_component + z_component);
    return kinetic_energy;
}


double compute_1e_V_ce_energy(MultiResolutionAnalysis<3> &MRA, FunctionTree<3> &Phi, FunctionTree<3> &V){
    FunctionTree<3> V_Phi(MRA);
    mrcpp::multiply(building_precision, V_Phi, 1.0, V, Phi);

    double V_ce = mrcpp::dot(Phi, V_Phi);
    return V_ce;
}


double compute_V_ee_energy_He(MultiResolutionAnalysis<3> &MRA, FunctionTree<3> &Phi){
    FunctionTree<3> J_Phi(MRA);
    FunctionTree<3> Rho(MRA);

    // Define density:
    mrcpp::multiply(building_precision, Rho, 1.0, Phi, Phi);

    // Build the convolution operator:
    mrcpp::PoissonOperator Poisson(MRA, building_precision);
    
    // Apply the convolution operator to get J|\Phi>:
    mrcpp::apply(building_precision, J_Phi, Poisson, Rho, -1, false);
    // J|\Phi> = \int d^3r' \frac{1}{|r-r'|} \rho(r') \Phi(r) --> This is missing, but we're integrating in the next step with respect to the density Rho.

    double V_ee = dot(Rho, J_Phi);

    return V_ee;
}



void compute_He_energy(MultiResolutionAnalysis<3> &MRA, FunctionTree<3> &Phi, FunctionTree<3> &V, double &mu, double &E){
    double T = compute_1e_Kinetic_energy(MRA, Phi);
    double V_ce = compute_1e_V_ce_energy(MRA, Phi, V);
    double V_ee = compute_V_ee_energy_He(MRA, Phi);
    double V_t = 2*V_ce + V_ee;
    double Exp_Val_F = T + V_ce + V_ee;

    std::cout << "> T = " << 2*T << '\n';           // T expected --> +2.85
    std::cout << "> V_ce = " << 2* V_ce << '\n';    // Vce expect --> -6.75
    std::cout << "> V_ee = " << V_ee << '\n';       // Vee expect --> +1.05
    std::cout << "> V_t = " << V_t << '\n';    



    double E_one_el = T + V_ce;
    E = 2 * E_one_el + V_ee;
    std::cout << '\n' << "$ E = " << E << '\n' <<'\n'; // -> -2.85516

    mu = std::sqrt(-2*Exp_Val_F);
}


void apply_Helmholtz_He(double &mu, MultiResolutionAnalysis<3> &MRA, FunctionTree<3> &Phi_out, FunctionTree<3> &Phi_in, FunctionTree<3> &V){
    
    HelmholtzOperator Helm(MRA, mu, building_precision);
    PoissonOperator Poisson(MRA, building_precision);


    FunctionTree<3> V_Phi(MRA);
    FunctionTree<3> Rho(MRA);
    FunctionTree<3> J_Phi_convolution(MRA);
    FunctionTree<3> J_Phi(MRA);

    

    // Start with V*Phi
    mrcpp::multiply(building_precision, V_Phi, -2.0, V, Phi_in);

    // Then the J_phi
    mrcpp::multiply(building_precision, Rho, 1.0, Phi_in, Phi_in);
    mrcpp::apply(building_precision, J_Phi_convolution, Poisson, Rho, -1, false);
    mrcpp::multiply(building_precision, J_Phi, -2.0, J_Phi_convolution, Phi_in);


    FunctionTree<3> To_convolute(MRA);
    mrcpp::add(building_precision, To_convolute, 1.0, V_Phi, 1.0, J_Phi);
    


    mrcpp::apply(building_precision, Phi_out, Helm, To_convolute, -1, false);
}
