#include "MRCPP/Gaussians"
#include "MRCPP/MWFunctions"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include <array>

const auto min_scale = -4;
const auto max_depth = 25;

const auto order = 7;
const auto prec = 1.0e-3;

const auto D = 3;

int main(int argc, char **argv) {
    auto timer = mrcpp::Timer();

    // Initialize printing
    auto printlevel = 0;
    mrcpp::Printer::init(printlevel);
    mrcpp::print::environment(0);
    mrcpp::print::header(0, "Multiplying MW functions");

    // Constructing world box
    auto corner = std::array<int, D>{-1, -1, -1};
    auto boxes = std::array<int, D>{2, 2, 2};
    auto world = mrcpp::BoundingBox<D>(min_scale, corner, boxes);

    // Constructing basis and MRA
    auto basis = mrcpp::InterpolatingBasis(order);
    auto MRA = mrcpp::MultiResolutionAnalysis<D>(world, basis, max_depth);

    // Setting up analytic Gaussians
    auto beta = 20.0;
    auto alpha = std::pow(beta / mrcpp::pi, 3.0 / 2.0);
    auto f_pos = mrcpp::Coord<D>{0.0, 0.0, 0.17};
    auto g_pos = mrcpp::Coord<D>{0.0, 0.0, -0.1};
    auto power = std::array<int, D>{0, 0, 0};
    auto f_func = mrcpp::GaussFunc<D>(beta, alpha, f_pos, power);
    auto g_func = mrcpp::GaussFunc<D>(beta, alpha, g_pos, power);

    // analytic product of two gaussians
    auto prod_pos =   mrcpp::Coord<D>{};
    auto beta_prod = beta*2;
    auto dist_2 = 0.0;
    for (int i=0;i< D;i++){
        dist_2 += (f_pos[i]-g_pos[i])*(f_pos[i]-g_pos[i]);//dot product of coordinate diff
        prod_pos[i] = 0.5*(f_pos[i]+g_pos[i]);
    }
    auto alpha_prod = alpha*alpha*exp(-beta*0.5*dist_2);
    auto prod_func = mrcpp::GaussFunc<D>(beta_prod, alpha_prod, prod_pos, power);

    // Initialize MW functions
    mrcpp::FunctionTree<D> f_tree(MRA);
    mrcpp::FunctionTree<D> g_tree(MRA);
    mrcpp::FunctionTree<D> h_tree(MRA);
    mrcpp::FunctionTree<D> h2_tree(MRA);
    mrcpp::FunctionTree<D> prod_tree(MRA);

    // Projecting f and g
    mrcpp::project<D>(prec, f_tree, f_func);
    mrcpp::project<D>(prec, g_tree, g_func);
    mrcpp::project<D>(prec/1000, prod_tree, prod_func);

    // h = f*g
    mrcpp::multiply(prec, h_tree, 1.0, f_tree, g_tree);
    mrcpp::multiply(prec, h2_tree, 1.0, f_tree, g_tree, -1, true, true);

    auto integral = h_tree.integrate();
    auto sq_norm = h_tree.getSquareNorm();
    mrcpp::print::value(0, "Integral method 1", integral);
    mrcpp::print::value(0, "Square norm method 1", sq_norm);

    integral = h2_tree.integrate();
    sq_norm = h2_tree.getSquareNorm();
    mrcpp::print::value(0, "Integral method 2", integral);
    mrcpp::print::value(0, "Square norm method 2", sq_norm);

    h_tree.add(-1.0,prod_tree);
    sq_norm = h_tree.getSquareNorm();
    mrcpp::print::value(0, "Square norm error method 1" , sq_norm);
    h2_tree.add(-1.0,prod_tree);
    sq_norm = h2_tree.getSquareNorm();
    mrcpp::print::value(0, "Square norm error method 2" , sq_norm);


    mrcpp::print::footer(0, timer, 2);

    return 0;
}
