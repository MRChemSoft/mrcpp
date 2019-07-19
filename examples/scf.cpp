#include "MRCPP/MWFunctions"
#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include <memory>

const auto min_scale = -4;
const auto max_depth = 25;

const auto order = 5;
const auto prec = 1.0e-3;
const auto D = 3; // Dimensions

using namespace mrcpp;

void setupNuclearPotential(double Z, FunctionTree<D> &V) {
    auto timer = Timer();
    auto oldlevel = Printer::setPrintLevel(10);
    print::header(0, "Projecting nuclear potential");

    // Smoothing parameter
    auto c = 0.00435 * prec / std::pow(Z, 5);
    auto u = [](double r) -> double {
        return std::erf(r) / r +
               1.0 / (3.0 * std::sqrt(mrcpp::pi)) * (std::exp(-r * r) + 16.0 * std::exp(-4.0 * r * r));
    };
    auto f = [u, c, Z](const Coord<3> &r) -> double {
        auto x = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
        return -1.0 * Z * u(x / c) / c;
    };

    // Projecting function
    project<D>(prec, V, f);

    print::footer(0, timer, 2);
    Printer::setPrintLevel(oldlevel);
}

void setupInitialGuess(FunctionTree<D> &phi) {
    auto timer = Timer();
    auto oldlevel = Printer::setPrintLevel(10);
    print::header(0, "Projecting initial guess");

    auto f = [](const Coord<D> &r) -> double {
        auto x = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
        return 1.0 * std::exp(-1.0 * x * x);
    };

    // Projecting and normalizing function
    project<D>(prec, phi, f);
    phi.normalize();

    print::footer(0, timer, 2);
    Printer::setPrintLevel(oldlevel);
}

int main(int argc, char **argv) {
    auto timer = Timer();

    // Initialize printing
    auto printlevel = 0;
    Printer::init(printlevel);
    print::environment(0);

    // Constructing world box
    auto min_scale = -4;
    auto corner = std::array<int, D>{-1, -1, -1};
    auto boxes = std::array<int, D>{2, 2, 2};
    auto world = BoundingBox<D>(min_scale, corner, boxes);
    // Constructing basis and MRA
    auto basis = InterpolatingBasis(order);
    auto MRA = MultiResolutionAnalysis<D>(world, basis, max_depth);

    // Nuclear potential
    auto Z = 1.0;
    FunctionTree<D> V(MRA);
    setupNuclearPotential(Z, V);

    // Wave function
    auto phi_n = std::make_shared<FunctionTree<D>>(MRA);
    auto phi_np1 = std::make_shared<FunctionTree<D>>(MRA);
    setupInitialGuess(*phi_n);

    print::header(0, "Running SCF");
    printout(0, " Iter");
    printout(0, "      E_np1          dE_n   ");
    printout(0, "   ||phi_np1||   ||dPhi_n||" << std::endl);
    print::separator(0, '-');

    // Orbtial energies
    auto epsilon_n = -0.5;
    auto epsilon_np1 = 0.0;
    auto d_epsilon_n = 0.0;

    auto iter = 1;
    auto error = 1.0;
    auto scf_t = std::vector<Timer>();
    while (error > 10 * prec) {
        Timer cycle_t;

        // Initialize Helmholtz operator
        if (epsilon_n > 0.0) epsilon_n *= -1.0;
        auto mu_n = std::sqrt(-2.0 * epsilon_n);
        HelmholtzOperator H(MRA, mu_n, prec);

        // Compute Helmholtz argument V*phi
        FunctionTree<D> Vphi(MRA);
        copy_grid(Vphi, *phi_n);                 // Copy grid from orbital
        multiply(prec, Vphi, 1.0, V, *phi_n, 1); // Relax grid max one level

        // Apply Helmholtz operator phi^n+1 = H[V*phi^n]
        apply(prec, *phi_np1, H, Vphi);
        phi_np1->rescale(-1.0 / (2.0 * mrcpp::pi));

        // Compute orbital residual
        FunctionTree<D> d_phi_n(MRA);
        copy_grid(d_phi_n, *phi_np1);                    // Copy grid from phi_np1
        add(-1.0, d_phi_n, 1.0, *phi_np1, -1.0, *phi_n); // No grid relaxation
        error = std::sqrt(d_phi_n.getSquareNorm());

        // Compute energy update <Vphi|d_phi>/||phi||
        d_epsilon_n = dot(Vphi, d_phi_n) / phi_np1->getSquareNorm();
        epsilon_np1 = epsilon_n + d_epsilon_n;

        int oldprec = Printer::getPrecision();
        printout(0, std::setw(3) << iter);
        Printer::setPrecision(10);
        printout(0, std::setw(19) << epsilon_np1);
        Printer::setPrecision(1);
        printout(0, std::setw(9) << d_epsilon_n);
        Printer::setPrecision(10);
        printout(0, std::setw(19) << phi_np1->getSquareNorm());
        Printer::setPrecision(1);
        printout(0, std::setw(9) << error);
        Printer::setPrecision(oldprec);
        printout(0, std::endl);

        // Prepare for next iteration
        epsilon_n = epsilon_np1;
        phi_n = phi_np1;
        phi_np1 = std::make_shared<FunctionTree<D>>(MRA);
        phi_n->normalize();

        cycle_t.stop();
        scf_t.push_back(cycle_t);
        iter++;
    }

    print::separator(0, '=', 2);

    print::header(0, "SCF timings");
    for (auto i = 0; i < scf_t.size(); i++) {
        std::stringstream o_cycle;
        o_cycle << "Time cycle " << std::setw(3) << i + 1;
        print::time(0, o_cycle.str(), scf_t[i]);
    }
    print::separator(0, '=', 2);
    print::header(0, "Final energy");
    print::value(0, "Eigenvalue", epsilon_n, "(au)");
    print::separator(0, '=', 2);

    return 0;
}
