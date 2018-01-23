#include "MRCPP/MWFunctions"
#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

const int min_scale = -4;
const int max_depth = 25;

const int order = 7;
const double prec = 1.0e-5;

int main(int argc, char **argv) {
    mrcpp::Timer timer;

    // Initialize printing
    int printlevel = 0;
    mrcpp::Printer::init(printlevel);
    mrcpp::Printer::printEnvironment();
    mrcpp::Printer::printHeader(0, "Applying derivative operator");

    // Constructing world box
    int corner[1] = {-1 };
    int boxes[1]  = { 2 };
    mrcpp::BoundingBox<1> world(min_scale, corner, boxes);

    // Constructing basis and MRA
    mrcpp::InterpolatingBasis basis(order);
    mrcpp::MultiResolutionAnalysis<1> MRA(world, basis, max_depth);

    // Setting up analytic functions
    const double alpha = 3.0;
    const double r_0 = mrcpp::pi - 3.0;
    auto f = [alpha, r_0] (const double *r) -> double {
        double R = fabs(r[0] - r_0);
        return exp(-alpha*R);
    };
    auto df = [alpha, r_0] (const double *r) -> double {
        double R = fabs(r[0] - r_0);
        double sign = 1.0;
        if (r[0] > r_0) sign = -1.0;
        return sign*alpha*exp(-alpha*R);
    };

    // Initializing MW functions and operators
    mrcpp::ABGVOperator<1> D(MRA, 0.0, 0.0);
    mrcpp::FunctionTree<1> f_tree(MRA);
    mrcpp::FunctionTree<1> df_tree(MRA);
    mrcpp::FunctionTree<1> dg_tree(MRA);
    mrcpp::FunctionTree<1> err_tree(MRA);

    // Projecting functions
    mrcpp::project(prec, f_tree, f);
    mrcpp::project(prec, df_tree, df);

    // Applying derivative operator
    mrcpp::apply(dg_tree, D, f_tree, 0);

    // Computing error
    mrcpp::add(-1.0, err_tree, 1.0, df_tree, -1.0, dg_tree);

    double f_int = f_tree.integrate();
    double f_norm = sqrt(f_tree.getSquareNorm());
    double df_int = df_tree.integrate();
    double df_norm = sqrt(df_tree.getSquareNorm());
    double dg_int = dg_tree.integrate();
    double dg_norm = sqrt(dg_tree.getSquareNorm());
    double abs_err = sqrt(err_tree.getSquareNorm());
    double rel_err = abs_err/df_norm;

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

