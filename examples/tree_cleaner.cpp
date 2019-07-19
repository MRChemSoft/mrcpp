#include "MRCPP/MWFunctions"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

const auto min_scale = -4;
const auto max_depth = 25;

const auto order = 7;
const auto prec = 1.0e-5;

const auto D = 3;
int main(int argc, char **argv) {
    auto timer = mrcpp::Timer();

    // Initialize printing
    auto printlevel = 0;
    mrcpp::Printer::init(printlevel);
    mrcpp::print::environment(0);
    mrcpp::print::header(0, "Testing TreeCleaner");

    // Constructing world box
    auto corner = std::array<int, D>{-1, -1, -1};
    auto boxes = std::array<int, D>{2, 2, 2};
    auto world = mrcpp::BoundingBox<D>(min_scale, corner, boxes);

    // Constructing basis and MRA
    auto basis = mrcpp::InterpolatingBasis(order);
    auto MRA = mrcpp::MultiResolutionAnalysis<D>(world, basis, max_depth);

    // Defining analytic function
    auto f = [](const mrcpp::Coord<D> &r) -> double {
        const auto beta = 100.0;
        const auto alpha = std::pow(beta / mrcpp::pi, 3.0 / 2.0);
        ;
        auto R = std::sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
        return alpha * std::exp(-beta * R * R);
    };

    // Initializing MW function
    mrcpp::FunctionTree<D> f_tree(MRA);

    auto iter = 0;
    auto n_nodes = 1;
    while (n_nodes > 0) {
        mrcpp::project<D>(-1.0, f_tree, f);         // Projecting on fixed grid
        n_nodes = mrcpp::refine_grid(f_tree, prec); // Refine grid
        mrcpp::clear_grid(f_tree);                  // Clear MW coefs
        printout(0, " iter " << std::setw(3) << iter++ << std::setw(45));
        printout(0, " n_nodes " << std::setw(5) << n_nodes << std::endl);
    }
    // Projecting on final converged grid
    mrcpp::project<D>(-1.0, f_tree, f);

    auto integral = f_tree.integrate();
    auto sq_norm = f_tree.getSquareNorm();
    mrcpp::print::separator(0, '-');
    mrcpp::print::value(0, "Integral", integral);
    mrcpp::print::value(0, "Square norm", sq_norm);
    mrcpp::print::footer(0, timer, 2);

    return 0;
}
