#include "MRCPP/Printer"
#include "MRCPP/Timer"
#include "MRCPP/Gaussians"
#include "MRCPP/Plotter"
#include "MRCPP/MWFunctions"
#include "MRCPP/MWOperators"

#include <cstdlib>

// This is your personal sandbox to test out ideas.
// See the examples/ directory for inspiration.
// Please do not commit your pilot code to git.

using namespace mrcpp;

// GLOBAL VARIABLES
int Z;

// ==================================================================================================================================

double energy(double precision, MultiResolutionAnalysis<3> &MRA, FunctionTree<3> &f, FunctionTree<3> &V){
    FunctionTree<3> f_V(MRA);
    FunctionTreeVector<3> Nabla_f; // THE GRADIENT WILL BE A VECTOR!!!!!!!
    double norm_sqrd = dot(f,f);
 
    mrcpp::ABGVOperator<3> D(MRA, 0.0, 0.0);
    
    // Compute the gradient of f
    Nabla_f = mrcpp::gradient(D,f);

    FunctionTree<3> Nabla_fx(MRA); 
    FunctionTree<3> Nabla_fy(MRA); 
    FunctionTree<3> Nabla_fz(MRA); 
    get_func(Nabla_f,0).deep_copy(&Nabla_fx);
    get_func(Nabla_f,1).deep_copy(&Nabla_fy);
    get_func(Nabla_f,2).deep_copy(&Nabla_fz);

    // Compute the product of the gradient and the function
    mrcpp::multiply(precision,f_V, 1.0, f, V);
    double potential_energy = dot(f,f_V) / norm_sqrd;

    // Compute the kinetic energy
    double nabla_f_norm_sqrd = dot(Nabla_fx,Nabla_fx) + dot(Nabla_fy,Nabla_fy) + dot(Nabla_fz,Nabla_fz);
    double kinetic_energy = 0.5000000000 * nabla_f_norm_sqrd/ norm_sqrd;
    
    double result = potential_energy + kinetic_energy;

    std::cout << "T = " << kinetic_energy << '\n';
    std::cout << "V = " << potential_energy << '\n' << '\n';

    std::cout << "  ########################### " << '\n';
    std::cout << "  | E = " << result << " |" << '\n';
    std::cout << "  ########################### " << '\n' << '\n';
    return result;
}

// ==================================================================================================================================

double energy_update(MultiResolutionAnalysis<3> &MRA, FunctionTree<3> &Tilde_phi, FunctionTree<3> &DeltaPhi, FunctionTree<3> &V){
    FunctionTree<3> V_DeltaPhi(MRA); // Allocate space for the product of the potential and the difference

    // norm = <Tilde_phi|Tilde_phi>
    double norm_sqrd = dot(Tilde_phi,Tilde_phi);
    
    
    // |V_DeltaPhi> = V * |DeltaPhi>
    mrcpp::multiply(1e-5,V_DeltaPhi, 1.0, V, DeltaPhi);
    
    // <Tilde_phi|V * DeltaPhi>
    double result = dot(Tilde_phi, V_DeltaPhi);
    // <Tilde_phi|V * DeltaPhi> / <Tilde_phi|Tilde_phi>
    return (result / norm_sqrd) ;
}


// ==================================================================================================================================



void apply_Helmholtz(double mu,double building_precision, MultiResolutionAnalysis<3> &MRA, FunctionTree<3> &out, FunctionTree<3> &GVPhi){
    HelmholtzOperator Helm(MRA, mu, building_precision);
    apply(building_precision, out, Helm, GVPhi, -1, false);
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
// ==================================================================================================================================
// ==================================================================================================================================

int main(int argc, char **argv) {
    Timer timer;

    // Leggi i parametri dal file
    auto parameters = readParameters("Input_Parameters.txt");

    // Assegna i valori ai parametri corrispondenti
    int order = static_cast<int>(parameters["order"]);
    int MaxLevel = static_cast<int>(parameters["MaxLevel"]);
    double building_precision = parameters["building_precision"];
    double epsilon = parameters["epsilon"];
    Z = static_cast<int>(parameters["Z"]);
    Z = 1;

    std::cout << "Parameters read from file: " << '\n';
    std::cout << "Order = " << order << '\n';
    std::cout << "MaxLevel = " << MaxLevel << '\n';
    std::cout << "Building precision = " << building_precision << '\n';
    std::cout << "Epsilon = " << epsilon << '\n';
    std::cout << "Z = " << Z << '\n';
    // Initialize printing
    int printlevel = 3;
    Printer::init(printlevel);
    print::environment(0);
    print::header(0, "MRCPP pilot code");

    // 1) We set the max order and max level of the MRA
    // order e MaxLevel sono già stati assegnati dai parametri letti

    // 2) Build the Box and MRA
    std::array<int,2> box = {-20,20};
    mrcpp::BoundingBox<3> bb(box);
    mrcpp::MultiResolutionAnalysis mra(bb, order, MaxLevel);
    // Also, we create the trees for each function we'll be dealing with
    mrcpp::FunctionTree<3> VPhi_tree(mra);  // -> This is the tree that will hold the function
    mrcpp::FunctionTree<3> GVPhi_tree(mra); // -> This is the tree that will hold the output of the convolution
    mrcpp::FunctionTree<3> Potential_tree(mra); // -> This is the tree that will hold the potential
    // Then the 2 auxiliary trees
    mrcpp::FunctionTree<3> Unormalized_difference_tree(mra);
    mrcpp::FunctionTree<3> Normalized_difference_tree(mra);
    mrcpp::FunctionTree<3> Phi_n_copy_tree(mra);

    // ==================================================================================================================================
    // ================= My function will be now a 1s H atom orbital in 3D, given by the following lambda function ======================
    // ==================================================================================================================================

    // 3) We define the operator G^\mu as the Helmltz Convolutional Operator:      [G^\mu]
    double mu = 1.0;
    //mrcpp::HelmholtzOperator Helm(mra, mu, building_precision);



    // 4) Define a trial function as a Gaussian 
    //   The function is defined as \Phi(r) = \alpha\exp(-\beta r^2)

    // Parameters:    
    double beta = 0.5;
    double alpha = (beta/3.14159265358979)*(beta/3.14159265358979)*(beta/3.14159265358979);
    alpha = std::sqrt(alpha);
    mrcpp::Coord<3> pos = {0,0,0};
    std::array<int,3> pow = {0,0,0};
    // Gaussian function
    mrcpp::GaussFunc<3> Phi_trial(beta, alpha, pos, pow);

    

    std::function<double(const Coord<3> &x)> Phi_t = [] (const mrcpp::Coord<3> &r) -> double {
        double beta = 3;
        double alpha = (beta/3.14159265358979)*(beta/3.14159265358979)*(beta/3.14159265358979);
        auto R = std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        return alpha * std::exp(-beta*R*R);
    };


    // 5) Let's DEFINE now the POTENTIAL function V(r) = -1/r
    /*
    std::function<double(const Coord<3> &x)> V = [](const Coord<3> &x) -> double {
        double r = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
        r = std::sqrt(r);
        double abs_r = std::abs(r);
        return -Z * 1.0 / abs_r;
    };
    */
   std::function<double(const Coord<3> &x)> V = [] (const mrcpp::Coord<3> &r) -> double {
        auto R = std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        return -Z/R;
    };



    // 6) PROJECT the function ON the TREE
    mrcpp::project<3>(building_precision, GVPhi_tree, Phi_trial); // I project the trial-function on the tree
    
//    mrcpp::project<3>(building_precision, GVPhi_tree, Phi_t); // I project the trial-function on the tree

    GVPhi_tree.normalize();
    mrcpp::project<3>(building_precision, Potential_tree, V); // I do the same with the potential


    // debug
    /*
    std::cout << "Potential_tree = " << '\n';
    std::cout << Potential_tree << '\n';


    std::cout << "Gauss_tree = " << '\n';
    std::cout << GVPhi_tree << '\n';
    */




    // ==================================================================================================================================
    // ================= Begin the SCF cycle ============================================================================================
    // ==================================================================================================================================

    // Parameters for the SCF
    double apply_precision = building_precision; // Precision of the convolution

    //int maxIter = -1; //          -> Max number of iterations for the convolution, if -1 it goes untill the precision is reached
    
    double norm_diff = 1;   //    -> Norm of the difference in the 2 consecutive iterations (initialized as 1 to begin the while loop)

    int num_cycle =  0; //         -> SCF cycle counter 

    double E = -1;  // Initialize the energy with a random number
    


    // CYCLE 0, TRIAL FUNCTION
    std::cout << "Cycle " << num_cycle << " done...  Norm of the difference = "<< norm_diff << '\n';
    ///E = energy(building_precision, mra, GVPhi_tree, Potential_tree);
    
    //std::cout << "Energy = " << E << '\n';
    mu = std::sqrt(-2.0 * E);
    std::cout << "mu = " << mu << '\n';
    std::cout << "************************************************************" << '\n';
    num_cycle++;

    

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

        // GVPhi_tree = \Tilde{\Phi^{n+1}}=  G^\mu (-2 * V \Phi^{n})
        apply_Helmholtz(mu,building_precision, mra, GVPhi_tree, VPhi_tree);
        GVPhi_tree.rescale(1.0 / (4 * 3.14159265358979));
        // Compute |\Tilde{\DeltaPhi}> = |\Tilde{GVPhi_tree}> - |Phi_n_copy_tree>
        mrcpp::add(apply_precision, Unormalized_difference_tree, +1.0, GVPhi_tree, -1.0, Phi_n_copy_tree);

        E = E + energy_update(mra, GVPhi_tree, Unormalized_difference_tree, Potential_tree);
        
    
        // Normalize the function: \Tilde{\Phi^{n+1}} --> \Phi^{n+1}
        GVPhi_tree.normalize();
    
        // Normalized_difference_tree = \Phi^{n+1} - \Phi^{n}
        mrcpp::add(apply_precision, Normalized_difference_tree, +1.0, GVPhi_tree, -1.0, Phi_n_copy_tree);
        // norm_diff = || \Phi^{n+1} - \Phi^{n} ||
        norm_diff = Normalized_difference_tree.getSquareNorm();
        norm_diff = std::sqrt(norm_diff);





        // Print the norm of the difference
        std::cout << "Cycle " << num_cycle << " done...  Norm of the difference = "<< norm_diff << '\n';
        
        //E = energy(building_precision, mra, GVPhi_tree, Potential_tree);
        std::cout << "  ########################### " << '\n';
        std::cout << "  | E = " << E << " |" << '\n';
        std::cout << "  ########################### " << '\n' << '\n';
        mu = std::sqrt(-2.0 * E);
        std::cout << "mu = " << mu << '\n';
        std::cout << '\n';
        std::cout << "************************************************************" << '\n';

        
        // Increment the cycle counter
        num_cycle++;
    }

    print::footer(0, timer, 2);
    return 0;
}