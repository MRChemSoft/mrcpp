#include "BoysFunction.h"
#include "trees/FunctionTree.h"
#include "core/InterpolatingBasis.h"
#include "utils/Printer.h"
#include "treebuilders/project.h"

namespace mrcpp {

BoysFunction::BoysFunction(int n, double p)
        : RepresentableFunction<1>(),
          order(n),
          prec(p),
          MRA(BoundingBox<1>(), InterpolatingBasis(13)) {
}

double BoysFunction::evalf(const double *r) const {
    int oldlevel = Printer::setPrintLevel(0);

    int n = this->order;
    double x = r[0];
    auto f = [x, n] (const double *t) -> double {
        double t_2 = t[0] * t[0];
        double xt_2 = x*t_2;
        double t_2n = 1.0;
        if (n > 0) {
            t_2n = pow(t_2, n);
        }
        return exp(-xt_2)*t_2n;
    };

    FunctionTree<1> tree(this->MRA);
    mrcpp::project(this->prec, tree, f);
    double result = tree.integrate();

    Printer::setPrintLevel(oldlevel);
    return result;
}

} // namespace mrcpp
