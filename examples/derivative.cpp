#include "MRCPP/MWFunctions"
#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

const int min_scale = -4;
const int max_scale = 20;

const int order = 7;
const double prec = 1.0e-5;

int main(int argc, char **argv) {
    Timer timer;

    // Initialize printing
    int printlevel = 0;
    Printer::init(printlevel);
    Printer::printEnvironment();
    Printer::printHeader(0, "Applying derivative operator");

    // Constructing world box
    int corner[1] = {-1 };
    int boxes[1]  = { 2 };
    BoundingBox<1> world(min_scale, corner, boxes);

    // Constructing basis and MRA
    InterpolatingBasis basis(order);
    MultiResolutionAnalysis<1> MRA(world, basis);

    // Setting up projector
    MWDerivative<1> apply(max_scale);
    MWProjector<1> project(prec, max_scale);
    MWAdder<1> add(-1.0, max_scale);

    // Setting up analytic functions
    const double alpha = 3.0;
    const double r_0 = pi - 3.0;
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
    ABGVOperator<1> D(MRA, 0.0, 0.0);
    FunctionTree<1> f_tree(MRA);
    FunctionTree<1> df_tree(MRA);
    FunctionTree<1> dg_tree(MRA);
    FunctionTree<1> err_tree(MRA);

    // Projecting functions
    project(f_tree, f);
    project(df_tree, df);

    // Applying derivative operator
    apply(dg_tree, D, f_tree, 0);

    // Computing error
    add(err_tree, 1.0, df_tree, -1.0, dg_tree);

    double f_int = f_tree.integrate();
    double f_norm = sqrt(f_tree.getSquareNorm());
    double df_int = df_tree.integrate();
    double df_norm = sqrt(df_tree.getSquareNorm());
    double dg_int = dg_tree.integrate();
    double dg_norm = sqrt(dg_tree.getSquareNorm());
    double abs_err = sqrt(err_tree.getSquareNorm());
    double rel_err = abs_err/df_norm;

    Printer::printSeparator(0, ' ');
    Printer::printDouble(0, "f_tree integral", f_int);
    Printer::printDouble(0, "f_tree norm", f_norm);
    Printer::printSeparator(0, ' ');
    Printer::printDouble(0, "df_tree integral", df_int);
    Printer::printDouble(0, "df_tree norm", df_norm);
    Printer::printSeparator(0, ' ');
    Printer::printDouble(0, "dg_tree integral", dg_int);
    Printer::printDouble(0, "dg_tree norm", dg_norm);
    Printer::printSeparator(0, ' ');
    Printer::printDouble(0, "absolute error", abs_err, 1);
    Printer::printDouble(0, "relative error", rel_err, 1);
    Printer::printSeparator(0, ' ');

    timer.stop();
    Printer::printFooter(0, timer, 2);

    return 0;
}
