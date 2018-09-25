#include "MRCPP/MWFunctions"
#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include <memory>

const auto min_scale = -4;
const auto max_depth = 25;

const auto order = 5;
const auto prec = 1.0e-2;

using namespace mrcpp;

void setupNuclearPotential(double Z, FunctionTree<3> &V) {
    Timer timer;
    auto oldlevel = Printer::setPrintLevel(10);
    Printer::printHeader(0, "Projecting nuclear potential");

    // Smoothing parameter
    auto c = 0.00435*prec/pow(Z, 5);
    auto u = [] (double r) -> double {
        return erf(r)/r + 1.0/(3.0*sqrt(mrcpp::pi))*(exp(-r*r) + 16.0*exp(-4.0*r*r));
    };
    auto f = [u, c, Z] (const double *r) -> double {
        auto x = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        return -1.0*Z*u(x/c)/c;
    };

    // Projecting function
    project(prec, V, f);

    timer.stop();
    Printer::printFooter(0, timer, 2);
    Printer::setPrintLevel(oldlevel);
}

void setupInitialGuess(FunctionTree<3> &phi) {
    Timer timer;
    auto oldlevel = Printer::setPrintLevel(10);
    Printer::printHeader(0, "Projecting initial guess");

    auto f = [] (const double *r) -> double {
        auto x = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        return 1.0*exp(-1.0*x*x);
    };

    // Projecting and normalizing function
    project(prec, phi, f);
    phi.normalize();

    timer.stop();
    Printer::printFooter(0, timer, 2);
    Printer::setPrintLevel(oldlevel);
}

int main(int argc, char **argv) {
    Timer timer;

    // Initialize printing
    auto printlevel = 0;
    Printer::init(printlevel);
    Printer::printEnvironment();

    // Constructing world box
    auto min_scale = -4;
    auto corner = std::vector<int>{-1, -1, -1};
    auto boxes = std::vector<int>{2, 2, 2};
    auto world = BoundingBox<3>(min_scale, corner, boxes);
    // Constructing basis and MRA
    auto basis = InterpolatingBasis(order);
    auto MRA = MultiResolutionAnalysis<3>(world, basis, max_depth);

    // Nuclear potential
    auto Z = 1.0;
    auto V = FunctionTree<3>(MRA);
    setupNuclearPotential(Z, V);

    // Wave function
    auto phi_n = std::make_shared<FunctionTree<3>>(MRA);
    auto phi_np1 = std::make_shared<FunctionTree<3>>(MRA);
    setupInitialGuess(*phi_n);

    Printer::printHeader(0, "Running SCF");
    printout(0, " Iter");
    printout(0, "      E_np1          dE_n   ");
    printout(0, "   ||phi_np1||   ||dPhi_n||" << std::endl);
    Printer::printSeparator(0, '-');

    // Orbtial energies
    auto epsilon_n = -0.5;
    auto epsilon_np1 = 0.0;
    auto d_epsilon_n = 0.0;

    auto iter = 1;
    auto error = 1.0;
    auto scf_t = std::vector<Timer>();
    while (error > 10*prec) {
        Timer cycle_t;

        // Initialize Helmholtz operator
        if (epsilon_n > 0.0) epsilon_n *= -1.0;
        auto mu_n = sqrt(-2.0*epsilon_n);
        auto H = HelmholtzOperator(MRA, mu_n, prec);

        // Compute Helmholtz argument V*phi
        auto Vphi = FunctionTree<3>(MRA);
        copy_grid(Vphi, *phi_n); // Copy grid from orbital
        multiply(prec, Vphi, 1.0, V, *phi_n, 1); // Relax grid max one level

        // Apply Helmholtz operator phi^n+1 = H[V*phi^n]
        apply(prec, *phi_np1, H, Vphi);
        phi_np1->rescale(-1.0/(2.0*mrcpp::pi));

        // Compute orbital residual
        auto d_phi_n = FunctionTree<3>(MRA);
        copy_grid(d_phi_n, *phi_np1); // Copy grid from phi_np1
        add(-1.0, d_phi_n, 1.0, *phi_np1, -1.0, *phi_n); // No grid relaxation
        error = sqrt(d_phi_n.getSquareNorm());

        // Compute energy update <Vphi|d_phi>/||phi||
        d_epsilon_n = dot(Vphi, d_phi_n)/phi_np1->getSquareNorm();
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

        // Prepare for next iteration
        epsilon_n = epsilon_np1;
        phi_n = phi_np1;
        phi_np1 = std::make_shared<FunctionTree<3>>(MRA);
        phi_n->normalize();

        cycle_t.stop();
        scf_t.push_back(cycle_t);
        iter++;
    }

    Printer::printSeparator(0, '=', 2);

    Printer::printHeader(0, "SCF timings");
    for (auto i = 0; i < scf_t.size(); i++) {
        Printer::printTree(0, "Time cycle", i+1, scf_t[i].getWallTime());
    }
    Printer::printSeparator(0, '=', 2);

    Printer::setPrecision(15);
    Printer::printHeader(0, "Final energy");
    Printer::printDouble(0, "Eigenvalue", epsilon_n);
    Printer::printSeparator(0, '=', 2);


    return 0;
}
