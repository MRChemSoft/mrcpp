#ifndef BOYSFUNCTION_H
#define BOYSFUNCTION_H

#include "RepresentableFunction.h"
#include "FunctionTree.h"
#include "MWProjector.h"
#include "InterpolatingBasis.h"
#include "constants.h"

class BoysFunction : public RepresentableFunction<1> {
public:
    BoysFunction(int n, double prec = 1.0e-10)
            : RepresentableFunction<1>(),
              Q(0),
              order(n) {
        int k = 13;
        InterpolatingBasis basis(k);
        BoundingBox<1> world;
        MultiResolutionAnalysis<1> MRA(world, basis);
        this->Q = new MWProjector<1>(MRA, prec);
    }
    virtual ~BoysFunction() {
        if (this->Q != 0) delete this->Q;
    }

    double evalf(const double *r) const {
        int oldlevel = TelePrompter::setPrintLevel(0);

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

        FunctionTree<1> *tree = (*this->Q)(f);
        double result = tree->integrate();
        delete tree;

        TelePrompter::setPrintLevel(oldlevel);
        return result;
    }

protected:
    MWProjector<1> *Q;
    const int order;
};

#endif // BOYSFUNCTION_H

