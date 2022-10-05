#include "MRCPP/MWFunctions"
#include "MRCPP/MWOperators"
#include "MRCPP/Gaussians"
#include "MRCPP/Printer"
#include "MRCPP/Plotter"
#include "MRCPP/Timer"
#include "MRCPP/utils/math_utils.h"

constexpr int D = 3;

extern int order;
extern int min_scale;
extern int max_depth;
extern double prec;
extern double bond_length;

mrcpp::MultiResolutionAnalysis<D> setup_mra() {
    // Constructing world box
    auto s_fac = std::array<double, 3>{2*bond_length, 2*bond_length, 2*bond_length};
    auto world = mrcpp::BoundingBox<3>(s_fac, true);

    // Constructing basis and MRA
    auto basis = mrcpp::InterpolatingBasis(order);
    return mrcpp::MultiResolutionAnalysis<D>(world, basis, max_depth);
}

void setup_plotter(mrcpp::Plotter<D> &plt) {
    double x = bond_length / 2.0;
    mrcpp::Coord<D> O{-3.0 * x, -x, -x};
    mrcpp::Coord<D> A{6.0 * x, 0.0, 0.0};
    plt.setOrigin(O);
    plt.setRange(A);
}

mrcpp::GaussExp<D> setup_monopole() {
    double c = bond_length;
    double x = bond_length / 2.0;
    mrcpp::GaussExp<D> f_exp;
    {
        // Setting up analytic Gaussian
        auto beta = 1.2;
        auto alpha = std::pow(beta / mrcpp::pi, 3.0 / 2.0);
        auto f_func = mrcpp::GaussFunc<D>(beta, alpha);
        f_func.setPos({c, c, c});
        f_exp.append(f_func);
    }
    return f_exp;
}

mrcpp::GaussExp<D> setup_dipole() {
    double c = bond_length;
    double x = bond_length / 2.0;
    mrcpp::GaussExp<D> f_exp;
    {
        // Setting up analytic Gaussian
        auto beta = 1.2;
        auto alpha = std::pow(beta / mrcpp::pi, 3.0 / 2.0);
        auto f_func = mrcpp::GaussFunc<D>(beta, alpha);
        f_func.setPos({c, c, c + x});
        f_exp.append(f_func);
    }
    {
        // Setting up analytic Gaussian
        auto beta = 0.6;
        auto alpha = -std::pow(beta / mrcpp::pi, 3.0 / 2.0);
        auto f_func = mrcpp::GaussFunc<D>(beta, alpha);
        f_func.setPos({c, c, c - x});
        f_exp.append(f_func);
    }
    return f_exp;
}

mrcpp::GaussExp<D> setup_quadrupole() {
    double c = bond_length;
    double x = bond_length / 2.0;
    mrcpp::GaussExp<D> f_exp;
    {
        // Setting up analytic Gaussian
        auto beta = 1.2;
        auto alpha = std::pow(beta / mrcpp::pi, 3.0 / 2.0);
        auto f_func = mrcpp::GaussFunc<D>(beta, alpha);
        f_func.setPos({c, c+x, c+x});
        f_exp.append(f_func);
        f_func.setPos({c, c-x, c-x});
        f_exp.append(f_func);
    }
    {
        // Setting up analytic Gaussian
        auto beta = 0.6;
        auto alpha = -std::pow(beta / mrcpp::pi, 3.0 / 2.0);
        auto f_func = mrcpp::GaussFunc<D>(beta, alpha);
        f_func.setPos({c, c-x, c+x});
        f_exp.append(f_func);
        f_func.setPos({c, c+x, c-x});
        f_exp.append(f_func);
    }
    return f_exp;
}

mrcpp::GaussExp<D> setup_ionic() {
    double c = bond_length;
    double x = bond_length / 2.0;
    mrcpp::GaussExp<D> f_exp;
    {
        // Setting up analytic Gaussian
        auto beta = 1.2;
        auto alpha = std::pow(beta / mrcpp::pi, 3.0 / 2.0);
        auto f_func = mrcpp::GaussFunc<D>(beta, alpha);
        f_func.setPos({c+x, c+x, c+x});
        f_exp.append(f_func);
        f_func.setPos({c-x, c-x, c+x});
        f_exp.append(f_func);
        f_func.setPos({c-x, c+x, c-x});
        f_exp.append(f_func);
        f_func.setPos({c+x, c-x, c-x});
        f_exp.append(f_func);
    }
    {
        // Setting up analytic Gaussian
        auto beta = 0.6;
        auto alpha = -std::pow(beta / mrcpp::pi, 3.0 / 2.0);
        auto f_func = mrcpp::GaussFunc<D>(beta, alpha);
        f_func.setPos({c-x, c+x, c+x});
        f_exp.append(f_func);
        f_func.setPos({c+x, c-x, c+x});
        f_exp.append(f_func);
        f_func.setPos({c+x, c+x, c-x});
        f_exp.append(f_func);
        f_func.setPos({c-x, c-x, c-x});
        f_exp.append(f_func);
    }
    return f_exp;
}

mrcpp::GaussExp<D> setup_covalent() {
    double c = bond_length;
    double x = bond_length / 2.0;
    mrcpp::GaussExp<D> f_exp;
    {
        // Setting up analytic Gaussian
        auto beta = 1.2;
        auto alpha = std::pow(beta / mrcpp::pi, 3.0 / 2.0);
        auto f_func = mrcpp::GaussFunc<D>(beta, alpha);
        f_func.setPos({c+x, c+x, c+x});
        f_exp.append(f_func);
        f_func.setPos({c-x, c+x, c+x});
        f_exp.append(f_func);
        f_func.setPos({c+x, c-x, c+x});
        f_exp.append(f_func);
        f_func.setPos({c+x, c+x, c-x});
        f_exp.append(f_func);
        f_func.setPos({c-x, c-x, c+x});
        f_exp.append(f_func);
        f_func.setPos({c-x, c+x, c-x});
        f_exp.append(f_func);
        f_func.setPos({c+x, c-x, c-x});
        f_exp.append(f_func);
        f_func.setPos({c-x, c-x, c-x});
        f_exp.append(f_func);
    }
    {
        // Setting up analytic Gaussian
        auto beta = 0.6;
        auto alpha = -std::pow(beta / mrcpp::pi, 3.0 / 2.0);
        auto f_func = mrcpp::GaussFunc<D>(beta, alpha);
        f_func.setPos({c+x, c+x, c+x});
        f_exp.append(f_func);
        f_func.setPos({c-x, c+x, c+x});
        f_exp.append(f_func);
        f_func.setPos({c+x, c-x, c+x});
        f_exp.append(f_func);
        f_func.setPos({c+x, c+x, c-x});
        f_exp.append(f_func);
        f_func.setPos({c-x, c-x, c+x});
        f_exp.append(f_func);
        f_func.setPos({c-x, c+x, c-x});
        f_exp.append(f_func);
        f_func.setPos({c+x, c-x, c-x});
        f_exp.append(f_func);
        f_func.setPos({c-x, c-x, c-x});
        f_exp.append(f_func);
    }
    return f_exp;
}

mrcpp::FunctionTree<D> * project_monopole() {
    // Setup environment
    auto MRA = setup_mra(0);
    auto *f_tree = new mrcpp::FunctionTree<D>(MRA);

    // Setup density
    auto s_fac = std::array<double, 3>{2.0 * bond_length, 2.0 * bond_length, 2.0 * bond_length};
    auto f_exp = setup_monopole();
    auto p_exp = f_exp.periodify(s_fac, 8.0);

    // Project function
    mrcpp::Timer t1;
    mrcpp::build_grid(*f_tree, p_exp);
    mrcpp::project(prec, *f_tree, p_exp);
    t1.stop();

    // Compute properties
    auto f_int = f_tree->integrate();
    auto f_norm = std::sqrt(f_tree->getSquareNorm());

    // Print properties
    println(0, std::setw(4) << "mono" <<
               std::setw(10) << std::setprecision(1) << t1.elapsed() <<
               std::setw(21) << " " <<
               std::setw(10) << " " <<
               std::setw(21) << std::setprecision(12) << f_norm <<
               std::setw(21) << std::setprecision(12) << f_int);
    return f_tree;
}

mrcpp::FunctionTree<D> * project_dipole() {
    // Setup environment
    auto MRA = setup_mra(0);
    auto *f_tree = new mrcpp::FunctionTree<D>(MRA);

    // Setup density
    auto s_fac = std::array<double, 3>{2.0 * bond_length, 2.0 * bond_length, 2.0 * bond_length};
    auto f_exp = setup_dipole();
    auto p_exp = f_exp.periodify(s_fac, 8.0);

    // Project function
    mrcpp::Timer t1;
    mrcpp::build_grid(*f_tree, p_exp);
    mrcpp::project(prec, *f_tree, p_exp);
    t1.stop();

    // Compute properties
    auto f_int = f_tree->integrate();
    auto f_norm = std::sqrt(f_tree->getSquareNorm());

    // Print properties
    println(0, std::setw(4) << "di" <<
               std::setw(10) << std::setprecision(1) << t1.elapsed() <<
               std::setw(21) << " " <<
               std::setw(10) << " " <<
               std::setw(21) << std::setprecision(12) << f_norm <<
               std::setw(21) << std::setprecision(12) << f_int);
    return f_tree;
}

mrcpp::FunctionTree<D> * project_quadrupole() {
    // Setup environment
    auto MRA = setup_mra(0);
    auto *f_tree = new mrcpp::FunctionTree<D>(MRA);

    // Setup density
    auto s_fac = std::array<double, 3>{2.0 * bond_length, 2.0 * bond_length, 2.0 * bond_length};
    auto f_exp = setup_quadrupole();
    auto p_exp = f_exp.periodify(s_fac, 8.0);

    // Project function
    mrcpp::Timer t1;
    mrcpp::build_grid(*f_tree, p_exp);
    mrcpp::project(prec, *f_tree, p_exp);
    t1.stop();

    // Compute properties
    auto f_int = f_tree->integrate();
    auto f_norm = std::sqrt(f_tree->getSquareNorm());

    // Print properties
    println(0, std::setw(4) << "quad" <<
               std::setw(10) << std::setprecision(1) << t1.elapsed() <<
               std::setw(21) << " " <<
               std::setw(10) << " " <<
               std::setw(21) << std::setprecision(12) << f_norm <<
               std::setw(21) << std::setprecision(12) << f_int);
    return f_tree;
}

mrcpp::FunctionTree<D> * project_ionic() {
    // Setup environment
    auto MRA = setup_mra(0);
    auto *f_tree = new mrcpp::FunctionTree<D>(MRA);

    // Setup density
    auto s_fac = std::array<double, 3>{2.0 * bond_length, 2.0 * bond_length, 2.0 * bond_length};
    auto f_exp = setup_ionic();
    auto p_exp = f_exp.periodify(s_fac, 8.0);

    // Project function
    mrcpp::Timer t1;
    mrcpp::build_grid(*f_tree, p_exp);
    mrcpp::project(prec, *f_tree, p_exp);
    t1.stop();

    // Compute properties
    auto f_int = f_tree->integrate();
    auto f_norm = std::sqrt(f_tree->getSquareNorm());

    // Print properties
    println(0, std::setw(4) << "ion" <<
               std::setw(10) << std::setprecision(1) << t1.elapsed() <<
               std::setw(21) << " " <<
               std::setw(10) << " " <<
               std::setw(21) << std::setprecision(12) << f_norm <<
               std::setw(21) << std::setprecision(12) << f_int);
    return f_tree;
}

mrcpp::FunctionTree<D> * project_covalent() {
    // Setup environment
    auto MRA = setup_mra(0);
    auto *f_tree = new mrcpp::FunctionTree<D>(MRA);

    // Setup density
    auto s_fac = std::array<double, 3>{2.0 * bond_length, 2.0 * bond_length, 2.0 * bond_length};
    auto f_exp = setup_covalent();
    auto p_exp = f_exp.periodify(s_fac, 8.0);

    // Project function
    mrcpp::Timer t1;
    mrcpp::build_grid(*f_tree, p_exp);
    mrcpp::project(prec, *f_tree, p_exp);
    t1.stop();

    // Compute properties
    auto f_int = f_tree->integrate();
    auto f_norm = std::sqrt(f_tree->getSquareNorm());

    // Print properties
    println(0, std::setw(4) << "cov" <<
               std::setw(10) << std::setprecision(1) << t1.elapsed() <<
               std::setw(21) << " " <<
               std::setw(10) << " " <<
               std::setw(21) << std::setprecision(12) << f_norm <<
               std::setw(21) << std::setprecision(12) << f_int);
    return f_tree;
}

void apply_operator(int n, std::vector<mrcpp::FunctionTree<D> *> &g_trees, mrcpp::FunctionTree<D> *f_tree, std::vector<double> &energies) {
    // Setup environment
    auto MRA = setup_mra(-n);
    mrcpp::PoissonOperator P(MRA, prec, -n, 5);
    auto *g_tree = new mrcpp::FunctionTree<D>(MRA);

    // Apply Poisson operator
    mrcpp::Timer t1;
    mrcpp::apply(prec, *g_tree, P, *f_tree);
    t1.stop();

    // Compute properties
    auto g_int = g_tree->integrate();
    auto g_norm = std::sqrt(g_tree->getSquareNorm());
    auto d_norm = 0.0;
    auto d_energy = 0.0;
    auto n_energy = mrcpp::dot(*g_tree, *f_tree);
    if (energies.size() > 0) d_energy = n_energy - energies.back();

    // Plot charge
    {
        mrcpp::FunctionTree<D> f_tmp(MRA);
        mrcpp::copy_grid(f_tmp, *f_tree);
        mrcpp::copy_func(f_tmp, *f_tree);
        mrcpp::refine_grid(f_tmp, 1);

        mrcpp::Plotter<D> plt;
        setup_plotter(plt);
        std::string f_name = "f_tree-" + std::to_string(n);
        plt.linePlot({10000}, f_tmp, f_name);
    }

    // Plot potential
    {
        mrcpp::FunctionTree<D> g_tmp(MRA);
        mrcpp::copy_grid(g_tmp, *g_tree);
        mrcpp::copy_func(g_tmp, *g_tree);
        mrcpp::refine_grid(g_tmp, 1);

        mrcpp::Plotter<D> plt;
        setup_plotter(plt);
        std::string g_name = "g_tree-" + std::to_string(n);
        plt.linePlot({10000}, g_tmp, g_name);
    }

    // Plot difference
    if (g_trees.size() > 0) {
        mrcpp::FunctionTree<D> d_tree(MRA);
        mrcpp::build_grid(d_tree, *g_tree);
        mrcpp::build_grid(d_tree, *g_trees.back());
        mrcpp::add(-1.0, d_tree, 1.0, *g_tree, -1.0, *g_trees.back());
        d_norm = d_tree.getSquareNorm();
        mrcpp::refine_grid(d_tree, 1);

        mrcpp::Plotter<D> plt;
        setup_plotter(plt);
        std::string d_name = "d_tree-" + std::to_string(n);
        plt.linePlot({10000}, d_tree, d_name);
    }

    // Print properties
    println(0, std::setw(4) << n <<
               std::setw(10) << std::setprecision(1) << t1.elapsed() <<
               std::setw(21) << std::setprecision(12) << n_energy <<
               std::setw(10) << std::setprecision(1) << d_energy <<
               std::setw(21) << std::setprecision(12) << g_norm <<
               std::setw(21) << std::setprecision(12) << d_norm);

    energies.push_back(n_energy);
    g_trees.push_back(g_tree);
}
