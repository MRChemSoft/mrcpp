#include "MRCPP/MWFunctions"
#include <MRCPP/MWOperators>
#include <MRCPP/Printer>
#include <MRCPP/Timer>
#include "operators/TimeEvolutionOperator.h"



std::complex<double> free_particle_analytical_solution(double x, double x0, double t, double sigma);

const auto min_scale = 0;
const auto max_depth = 25;

const auto order = 7;
const auto prec = 1.0e-7;

int main(int argc, char **argv)
{
    auto timer = mrcpp::Timer();

    // Initialize printing
    auto printlevel = 0;
    mrcpp::Printer::init(printlevel);
    mrcpp::print::environment(0);
    mrcpp::print::header(0, "Building operator");

    // Initialize world in the unit cube [0,1]
    auto basis = mrcpp::LegendreBasis(order);
    auto world = mrcpp::BoundingBox<1>(min_scale);
    auto MRA = mrcpp::MultiResolutionAnalysis<1>(world, basis, max_depth);

    mrcpp::print::footer(0, timer, 2);

    // Time moment:
    double t = 0.001;
    int finest_scale = 8;

    mrcpp::TimeEvolutionOperator<1> ReExp(MRA, prec, t, finest_scale, false);
    //Exp.initialize(4);  //impossible to initialise outside of a constructor!?
    println(0, ReExp.getComponent(0, 0));

    mrcpp::TimeEvolutionOperator<1> ImExp(MRA, prec, t, finest_scale, true);
    println(0, ImExp.getComponent(0, 0));

    //double x = 1.5;  // Example value for x
    //double result = std::erf(x);
    //std::cout << "erf(" << x << ") = " << result << std::endl;

    // Analytical solution
    double sigma = t;
    double x0 = 0.5;
    auto Re_f = [sigma, x0, t=0](const mrcpp::Coord<1> &r) -> double
    {
        return free_particle_analytical_solution(r[0], x0, t, sigma).real();
    };
    /*
    auto Im_f = [sigma, x0, t=0](const mrcpp::Coord<1> &r) -> double
    {
        return free_particle_analytical_solution(r[0], x0, t, sigma).imag();
    };
    */
//    std::cout << "t = " << t << std::endl;
    auto Re_g = [sigma, x0, t](const mrcpp::Coord<1> &r) -> double
    {
        return free_particle_analytical_solution(r[0], x0, t, sigma).real();
    };
    auto Im_g = [sigma, x0, t](const mrcpp::Coord<1> &r) -> double
    {
        return free_particle_analytical_solution(r[0], x0, t, sigma).imag();
    };
    // Test the functions
    /*
    for( double x = 0.0; x < 1.0; x += 0.1 )
        std::cout << "f(" << x << ") = " << Re_f({x}) << std::endl;
    for( double x = 0.0; x < 1.0; x += 0.1 )
        std::cout << "g(" << x << ") = " << Re_g({x}) << std::endl;
    */
    // Projecting functions
    mrcpp::FunctionTree<1> Re_f_tree(MRA);
    mrcpp::project<1>(prec, Re_f_tree, Re_f, -1);
//    auto integral = Re_f_tree.integrate();
    mrcpp::FunctionTree<1> Re_g_tree(MRA);
    mrcpp::project<1>(prec, Re_g_tree, Re_g, -1);
    mrcpp::FunctionTree<1> Im_g_tree(MRA);
    mrcpp::project<1>(prec, Im_g_tree, Im_g, -1);

    //std::cout << Re_f_tree;
    //std::cout << Re_g_tree;
    //std::cout << Im_g_tree;

    mrcpp::FunctionTree<1> Re_fout_tree(MRA);
    mrcpp::FunctionTree<1> Im_fout_tree(MRA);
    //mrcpp::FunctionTree<1> temp_fout_tree(MRA);
    mrcpp::apply(prec, Re_fout_tree, ReExp, Re_f_tree);
    //mrcpp::apply(prec, temp_fout_tree, ImExp, Im_f_tree);
    //Re_fout_tree -= temp_fout_tree;
    mrcpp::apply(prec, Im_fout_tree, ImExp, Re_f_tree);
    //... in case Im_f is not zero

    mrcpp::FunctionTree<1> Re_error(MRA);  // = Re_fout_tree - Re_g_tree
    mrcpp::FunctionTree<1> Im_error(MRA);  // = Im_fout_tree - Im_g_tree
    
    // Re_error = Re_fout_tree - Re_g_tree
    auto Re_f_vec = mrcpp::FunctionTreeVector<1>();
    Re_f_vec.push_back(std::make_tuple(1.0, &Re_fout_tree));
    Re_f_vec.push_back(std::make_tuple(-1.0, &Re_g_tree));
    mrcpp::add(prec, Re_error, Re_f_vec);

    auto Re_integral = Re_error.integrate();
    auto Re_sq_norm = Re_error.getSquareNorm();
    mrcpp::print::value(0, "Integral", Re_integral);
    mrcpp::print::value(0, "Square norm", Re_sq_norm);
    
    // Im_error = Im_fout_tree - Im_g_tree
    auto Im_f_vec = mrcpp::FunctionTreeVector<1>();
    Im_f_vec.push_back(std::make_tuple(1.0, &Im_fout_tree));
    Im_f_vec.push_back(std::make_tuple(-1.0, &Im_g_tree));
    mrcpp::add(prec, Im_error, Im_f_vec);

    auto Im_integral = Im_error.integrate();
    auto Im_sq_norm = Im_error.getSquareNorm();
    mrcpp::print::value(0, "Integral", Im_integral);
    mrcpp::print::value(0, "Square norm", Im_sq_norm);
    
    //std::cout << Re_error;
    //Re_error.crop(prec);
    //std::cout << Re_error;
    
    

    mrcpp::print::footer(0, timer, 2);
    return 0;
}


std::complex<double> free_particle_analytical_solution(double x, double x0, double t, double sigma)
{
    std::complex<double> i(0.0, 1.0);  // Imaginary unit
    auto denominator = 4 * t * i + sigma;
    std::complex<double> sqrt_denom = std::sqrt(denominator);
    std::complex<double> exponent = -((x - x0) * (x - x0)) / denominator;

    return std::sqrt(sigma) / sqrt_denom * std::exp(exponent);
}
