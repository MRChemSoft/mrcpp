#include "pbc.h"

int order = 5;
int min_scale = 0;
int max_depth = 25;
double prec = 1.0e-4;
double bond_length = 5.3;
std::string type = "ionic";

int main(int argc, char **argv) {
    auto timer = mrcpp::Timer();

    if (argc > 1) type = argv[1];
    if (argc > 2) prec = std::atof(argv[2]);
    if (argc > 3) order = std::atoi(argv[3]);
    if (argc > 4) min_scale = std::atoi(argv[4]);
    if (argc > 5) MSG_ABORT("Too many arguments!");

    auto printlevel = 0;
    mrcpp::Printer::init(printlevel);
    mrcpp::Printer::setWidth(90);
    mrcpp::print::environment(0);

    println(0, "type  : " << type);
    println(0, "prec  : " << prec);
    println(0, "order : " << order);
    println(0, "scale : " << -min_scale);
    mrcpp::print::separator(0, ' ', 1);

    mrcpp::FunctionTree<D> *f_tree = nullptr;
    std::vector<mrcpp::FunctionTree<D> *> g_trees;
    std::vector<double> energies;

    mrcpp::print::header(0, "Projecting charge density");
    println(0, std::setw(4) << "f" <<
               std::setw(10) << "t_proj" <<
               std::setw(21) << " " <<
               std::setw(10) << " " <<
               std::setw(21) << "sq_norm" <<
               std::setw(21) << "integral");
    mrcpp::print::separator(0, '-');
    if (type == "ionic") f_tree = project_ionic();
    if (type == "covalent") f_tree = project_covalent();
    if (type == "monopole") f_tree = project_monopole();
    if (type == "dipole") f_tree = project_dipole();
    if (type == "quadrupole") f_tree = project_quadrupole();
    mrcpp::print::footer(0, timer, 2);

    mrcpp::print::header(0, "Applying Poisson operator");
    println(0, std::setw(4) << "n" <<
               std::setw(10) << "t_apply" <<
               std::setw(21) << "energy" <<
               std::setw(10) << "update" <<
               std::setw(21) << "sq_norm" <<
               std::setw(21) << "diff_norm");
    mrcpp::print::separator(0, '-');
    for (int n = min_scale; n <= min_scale; n++) apply_operator(n, g_trees, f_tree, energies);
    mrcpp::print::footer(0, timer, 2);

    delete f_tree;
    for (int i = 0; i < g_trees.size(); i++) delete g_trees[i];

    return 0;
}
