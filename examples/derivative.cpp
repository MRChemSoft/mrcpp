#include "MRCPP/MWFunctions"
#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

const auto min_scale = -4;
const auto max_depth = 25;

const auto order = 7;
const auto prec = 1.0e-5;

const auto D = 1;

int main(int argc, char **argv) {
    auto timer = mrcpp::Timer();

    // Initialize printing
    auto printlevel = 0;
    mrcpp::Printer::init(printlevel);
    mrcpp::Printer::printEnvironment();
    mrcpp::Printer::printHeader(0, "Applying derivative operator");

    // Constructing world box
    auto corner = std::array<int, D>{-1};
    auto boxes = std::array<int, D>{2};
    auto world = mrcpp::BoundingBox<D>(min_scale, corner, boxes);

    // Constructing basis and MRA
    auto basis = mrcpp::InterpolatingBasis(order);
    auto MRA = mrcpp::MultiResolutionAnalysis<D>(world, basis, max_depth);

    // Setting up analytic functions
    const auto alpha = 3.0;
    const auto r_0 = mrcpp::pi - 3.0;
    auto f = [alpha, r_0](const mrcpp::Coord<D> &r) -> double {
        auto R = std::abs(r[0] - r_0);
        return std::exp(-alpha * R);
    };
    auto df = [alpha, r_0](const mrcpp::Coord<D> &r) -> double {
        auto R = std::abs(r[0] - r_0);
        auto sign = 1.0;
        if (r[0] > r_0) sign = -1.0;
        return sign * alpha * std::exp(-alpha * R);
    };

    // Initializing MW functions and operators
    mrcpp::ABGVOperator<D> D_00(MRA, 0.0, 0.0);
    mrcpp::FunctionTree<D> f_tree(MRA);
    mrcpp::FunctionTree<D> df_tree(MRA);
    mrcpp::FunctionTree<D> dg_tree(MRA);
    mrcpp::FunctionTree<D> err_tree(MRA);

    // Projecting functions
    mrcpp::project<D>(prec, f_tree, f);
    mrcpp::project<D>(prec, df_tree, df);

    // Applying derivative operator
    mrcpp::apply(dg_tree, D_00, f_tree, 0);

    // Computing error
    mrcpp::add(-1.0, err_tree, 1.0, df_tree, -1.0, dg_tree);

    auto f_int = f_tree.integrate();
    auto f_norm = std::sqrt(f_tree.getSquareNorm());
    auto df_int = df_tree.integrate();
    auto df_norm = std::sqrt(df_tree.getSquareNorm());
    auto dg_int = dg_tree.integrate();
    auto dg_norm = std::sqrt(dg_tree.getSquareNorm());
    auto abs_err = std::sqrt(err_tree.getSquareNorm());
    auto rel_err = abs_err / df_norm;

    mrcpp::Printer::printSeparator(0, ' ');
    mrcpp::Printer::printDouble(0, "f_tree integral", f_int);
    mrcpp::Printer::printDouble(0, "f_tree norm", f_norm);
    mrcpp::Printer::printSeparator(0, ' ');
    mrcpp::Printer::printDouble(0, "df_tree integral", df_int);
    mrcpp::Printer::printDouble(0, "df_tree norm", df_norm);
    mrcpp::Printer::printSeparator(0, ' ');
    mrcpp::Printer::printDouble(0, "dg_tree integral", dg_int);
    mrcpp::Printer::printDouble(0, "dg_tree norm", dg_norm);
    mrcpp::Printer::printSeparator(0, ' ');
    mrcpp::Printer::printDouble(0, "absolute error", abs_err, 1);
    mrcpp::Printer::printDouble(0, "relative error", rel_err, 1);
    mrcpp::Printer::printSeparator(0, ' ');

    timer.stop();
    mrcpp::Printer::printFooter(0, timer, 2);

    return 0;
}
