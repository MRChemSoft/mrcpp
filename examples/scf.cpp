#include "MRCPP/MWFunctions"
#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

const int min_scale = -4;
const int max_scale = 20;

const int order = 7;
const double prec = 1.0e-5;

void setupNuclearPotential(double Z, FunctionTree<3> &V) {
    Timer timer;
    int oldlevel = Printer::setPrintLevel(10);
    Printer::printHeader(0, "Projecting nuclear potential");

    // Setting up projector
    MWProjector<3> project(prec, max_scale);

    // Smoothing parameter
    double c = 0.00435*prec/pow(Z, 5);
    auto u = [] (double r) -> double {
        return erf(r)/r + 1.0/(3.0*sqrt(pi))*(exp(-r*r) + 16.0*exp(-4.0*r*r));
    };
    auto f = [u, c, Z] (const double *r) -> double {
        double x = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        return -1.0*Z*u(x/c)/c;
    };

    // Projecting function
    project(V, f);

    timer.stop();
    Printer::printFooter(0, timer, 2);
    Printer::setPrintLevel(oldlevel);
}

void setupInitialGuess(FunctionTree<3> &phi) {
    Timer timer;
    int oldlevel = Printer::setPrintLevel(10);
    Printer::printHeader(0, "Projecting initial guess");

    // Setting up projector
    MWProjector<3> project(prec, max_scale);

    auto f = [] (const double *r) -> double {
        double x = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        return 1.0*exp(-1.0*x*x);
    };

    // Projecting and normalizing function
    project(phi, f);
    phi.normalize();

    timer.stop();
    Printer::printFooter(0, timer, 2);
    Printer::setPrintLevel(oldlevel);
}

int main(int argc, char **argv) {
    Timer timer;

    // Initialize printing
    int printlevel = 0;
    Printer::init(printlevel);
    Printer::printEnvironment();

    // Constructing world box
    int min_scale = -4;
    int corner[3] = {-1,-1,-1};
    int boxes[3]  = { 2, 2, 2};
    NodeIndex<3> idx(min_scale, corner);
    BoundingBox<3> world(idx, boxes);

    // Constructing basis and MRA
    InterpolatingBasis basis(order);
    MultiResolutionAnalysis<3> MRA(world, basis);

    MWAdder<3> add(-1.0, max_scale);
    MWMultiplier<3> mult(prec, max_scale);
    MWConvolution<3> apply(prec, max_scale);
    GridGenerator<3> grid(max_scale);

    // Nuclear potential
    double Z = 1.0;
    FunctionTree<3> V(MRA);
    setupNuclearPotential(Z, V);

    // Wave function
    FunctionTree<3> *phi_n = new FunctionTree<3>(MRA);
    FunctionTree<3> *phi_np1 = 0;
    setupInitialGuess(*phi_n);

    Printer::printHeader(0, "Running SCF");
    printout(0, " Iter");
    printout(0, "      E_np1          dE_n   ");
    printout(0, "   ||phi_np1||   ||dPhi_n||" << std::endl);
    Printer::printSeparator(0, '-');

    // Orbtial energies
    double epsilon_n = -0.5;
    double epsilon_np1 = 0.0;
    double d_epsilon_n = 0.0;

    int iter = 1;
    double error = 1.0;
    std::vector<Timer> scf_t;
    while (error > 10*prec) {
        Timer cycle_t;

        // Initialize Helmholtz operator
        if (epsilon_n > 0.0) epsilon_n *= -1.0;
        double mu_n = sqrt(-2.0*epsilon_n);
        HelmholtzOperator H(MRA, mu_n, prec);

        // Compute Helmholtz argument V*phi
        FunctionTree<3> Vphi(MRA);
        grid(Vphi, *phi_n); // Copy grid from orbital
        mult(Vphi, 1.0, V, *phi_n, 1); // Relax grid max one level

        // Apply Helmholtz operator phi^n+1 = H[V*phi^n]
        phi_np1 = new FunctionTree<3>(MRA);
        apply(*phi_np1, H, Vphi);
        *phi_np1 *= -1.0/(2.0*pi);

        // Compute orbital residual
        FunctionTree<3> d_phi_n(MRA);
        grid(d_phi_n, *phi_np1); // Copy grid from phi_np1
        add(d_phi_n, 1.0, *phi_np1, -1.0, *phi_n); // No grid relaxation
        error = sqrt(d_phi_n.getSquareNorm());

        // Compute energy update <Vphi|d_phi>/||phi||
        d_epsilon_n = Vphi.dot(d_phi_n)/phi_np1->getSquareNorm();
        epsilon_np1 = epsilon_n + d_epsilon_n;

        printout(0, std::setw(3) << iter);
        Printer::setPrecision(10);
        printout(0, std::setw(19) << epsilon_np1);
        Printer::setPrecision(1);
        printout(0, std::setw(9) << d_epsilon_n);
        Printer::setPrecision(10);
        printout(0, std::setw(19) << phi_np1->getSquareNorm());
        Printer::setPrecision(1);
        printout(0, std::setw(9) << error);
        Printer::setPrecision(15);
        printout(0, std::endl);

        delete phi_n;

        // Prepare for next iteration
        epsilon_n = epsilon_np1;
        phi_n = phi_np1;
        phi_n->normalize();

        cycle_t.stop();
        scf_t.push_back(cycle_t);
        iter++;
    }
    Printer::printSeparator(0, '=', 2);

    Printer::printHeader(0, "SCF timings");
    for (int i = 0; i < scf_t.size(); i++) {
        println(0, " Time cycle " << std::setw(5) << i+1 << "  " << scf_t[i]);
    }
    Printer::printSeparator(0, '=', 2);

    Printer::setPrecision(15);
    Printer::printHeader(0, "Final energy");
    println(0, " Eigenvalue        " << std::setw(40) << epsilon_n);
    Printer::printSeparator(0, '=', 2);

    delete phi_n;

    return 0;
}

