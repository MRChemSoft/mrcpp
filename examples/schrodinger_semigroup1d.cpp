#include "MRCPP/MWFunctions"
#include "MRCPP/Plotter"
#include "functions/special_functions.h"
#include "operators/TimeEvolutionOperator.h"
#include "treebuilders/complex_apply.h"
#include <MRCPP/MWOperators>
#include <MRCPP/Printer>
#include <MRCPP/Timer>

const auto min_scale = 0;
const auto max_depth = 25;

const auto order = 4;
const auto prec = 1.0e-7;

int finest_scale = 10; // for time evolution operator construction (not recommended to use more than 10)
int max_Jpower = 20;   // the amount of J integrals to be used in construction (20 should be enough)

// Time moments:
double t1 = 0.001;        // initial time moment (not recommended to use more than 0.001)
double delta_t = 0.001;   // time step (not recommended to use less than 0.001)
double t2 = delta_t + t1; // final time moment

/**
 * @brief Exploring free-particle time evolution.
 * @details We check the time propagator.
 *
 * The time evolution equation is given by:
 * \f[
 *   g(x) = \exp \left( i \delta t \partial_x^2 \right) f(x)
 * \f]
 *
 * where \f$f(x) = \psi(x, t_1)\f$,
 * \f$g(x) = \psi(x, t_2)\f$,
 * \f$\delta t = t_2 - t_1\f$
 * and
 * \f[
 *   \psi(x, t) = \sqrt{\frac{\sigma}{4it + \sigma}} e^{-\frac{(x - x_0)^2}{4it + \sigma}}
 *   .
 * \f]
 *
 */
int main(int argc, char **argv) {
    auto timer = mrcpp::Timer();

    // Initialize printing
    auto printlevel = 0;
    mrcpp::Printer::init(printlevel);
    mrcpp::print::environment(0);

    // Initialize world in the unit cube [0,1]
    auto basis = mrcpp::LegendreBasis(order);
    auto world = mrcpp::BoundingBox<1>(min_scale);
    auto MRA = mrcpp::MultiResolutionAnalysis<1>(world, basis, max_depth);

    mrcpp::print::header(0, "Building operator");
    mrcpp::print::footer(0, timer, 2);

    // Time evolution operatror Exp(delta_t)
    mrcpp::TimeEvolutionOperator<1> ReExp(MRA, prec, delta_t, finest_scale, false);
    println(0, ReExp.getComponent(0, 0));
    mrcpp::TimeEvolutionOperator<1> ImExp(MRA, prec, delta_t, finest_scale, true);
    println(0, ImExp.getComponent(0, 0));

    mrcpp::print::header(0, "Preparing analytical solution");
    mrcpp::print::footer(0, timer, 2);

    // Analytical solution parameters for psi(x, t)
    double sigma = 0.001;
    double x0 = 0.5;

    // Functions f(x) = psi(x, t1) and g(x) = psi(x, t2)
    auto Re_f = [sigma, x0, t = t1](const mrcpp::Coord<1> &r) -> double { return mrcpp::free_particle_analytical_solution(r[0], x0, t, sigma).real(); };
    auto Im_f = [sigma, x0, t = t1](const mrcpp::Coord<1> &r) -> double { return mrcpp::free_particle_analytical_solution(r[0], x0, t, sigma).imag(); };
    auto Re_g = [sigma, x0, t = t2](const mrcpp::Coord<1> &r) -> double { return mrcpp::free_particle_analytical_solution(r[0], x0, t, sigma).real(); };
    auto Im_g = [sigma, x0, t = t2](const mrcpp::Coord<1> &r) -> double { return mrcpp::free_particle_analytical_solution(r[0], x0, t, sigma).imag(); };

    // Projecting functions
    mrcpp::FunctionTree<1> Re_f_tree(MRA);
    mrcpp::project<1, double>(prec, Re_f_tree, Re_f);
    mrcpp::FunctionTree<1> Im_f_tree(MRA);
    mrcpp::project<1, double>(prec, Im_f_tree, Im_f);
    mrcpp::FunctionTree<1> Re_g_tree(MRA);
    mrcpp::project<1, double>(prec, Re_g_tree, Re_g);
    mrcpp::FunctionTree<1> Im_g_tree(MRA);
    mrcpp::project<1, double>(prec, Im_g_tree, Im_g);

    // Output function trees
    mrcpp::FunctionTree<1> Re_fout_tree(MRA);
    mrcpp::FunctionTree<1> Im_fout_tree(MRA);

    // Complex objects for use in apply()
    mrcpp::ComplexObject<mrcpp::ConvolutionOperator<1>> E(ReExp, ImExp);
    mrcpp::ComplexObject<mrcpp::FunctionTree<1>> input(Re_f_tree, Im_f_tree);
    mrcpp::ComplexObject<mrcpp::FunctionTree<1>> output(Re_fout_tree, Im_fout_tree);

    mrcpp::print::header(0, "Applying operator");
    mrcpp::print::footer(0, timer, 2);

    // Apply operator Exp(delta_t) f(x)
    mrcpp::apply(prec, output, E, input);

    mrcpp::print::header(0, "Checking the result on analytical solution");
    mrcpp::print::footer(0, timer, 2);

    // Check g(x) = Exp(delta_t) f(x)
    mrcpp::FunctionTree<1> Re_error(MRA); // = Re_fout_tree - Re_g_tree
    mrcpp::FunctionTree<1> Im_error(MRA); // = Im_fout_tree - Im_g_tree

    // Re_error = Re_fout_tree - Re_g_tree
    add(prec, Re_error, 1.0, Re_fout_tree, -1.0, Re_g_tree);
    auto Re_integral = Re_error.integrate();
    auto Re_sq_norm = Re_error.getSquareNorm();
    mrcpp::print::value(0, "Integral of    Re(Exp(delta_t) f(x) - g(x)) =", Re_integral);
    mrcpp::print::value(0, "Square norm of Re(Exp(delta_t) f(x) - g(x)) =", Re_sq_norm);

    // Im_error = Im_fout_tree - Im_g_tree
    add(prec, Im_error, 1.0, Im_fout_tree, -1.0, Im_g_tree);
    auto Im_integral = Im_error.integrate();
    auto Im_sq_norm = Im_error.getSquareNorm();
    mrcpp::print::value(0, "Integral of    Im(Exp(delta_t) f(x) - g(x)) =", Im_integral);
    mrcpp::print::value(0, "Square norm of Im(Exp(delta_t) f(x) - g(x)) =", Im_sq_norm);

    mrcpp::print::header(0, "Saving plots to files");
    mrcpp::print::footer(0, timer, 2);

    // Set plotting parameters
    int nPts = 1000;
    mrcpp::Coord<1> o{0.0};
    mrcpp::Coord<1> a{1.0};
    mrcpp::Plotter<1> plot(o);
    plot.setRange(a);

    plot.linePlot({nPts}, Re_error, "Re_error");   // Write to file Re_error.line
    plot.linePlot({nPts}, Im_error, "Im_error");   // Write to file Im_error.line
    plot.linePlot({nPts}, Re_f_tree, "Re_f_tree"); // Write to file Re_f_tree.line
    plot.linePlot({nPts}, Im_f_tree, "Im_f_tree"); // Write to file Im_f_tree.line
    plot.linePlot({nPts}, Re_g_tree, "Re_g_tree"); // Write to file Re_g_tree.line
    plot.linePlot({nPts}, Im_g_tree, "Im_g_tree"); // Write to file Im_g_tree.line

    mrcpp::print::footer(0, timer, 2);
    return 0;
}
