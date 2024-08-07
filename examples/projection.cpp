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

    // Projecting function
    mrcpp::FunctionTree<D> f_tree(MRA);
    mrcpp::project<D, double>(prec, f_tree, f, -1);
    auto integral = f_tree.integrate();

    mrcpp::print::header(0, "Projecting analytic function");
    mrcpp::print::tree(-1, "Projected function", f_tree, timer);
    mrcpp::print::value(0, "Integrated function", integral, "(au)");
    mrcpp::print::footer(0, timer, 2);

    return 0;
}
