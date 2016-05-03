#ifndef BOYSFUNCTION_H
#define BOYSFUNCTION_H

#include "RepresentableFunction.h"
#include "FunctionTree.h"
#include "MWProjector.h"
#include "InterpolatingBasis.h"
#include "constants.h"

class BoysIntegrand : public RepresentableFunction<1> {
public:
    BoysIntegrand(int n, double x)
            : order(n),
              pos(x) {
        if (this->order < 0) MSG_FATAL("BoysFunction of negative order");
    }
    ~BoysIntegrand() { }

    double evalf(const double *r) const {
        double t_2 = (*r) * (*r);
        double x = this->pos;
        double xt_2 = x*t_2;
        double t_2n = 1.0;
        if (this->order > 0) {
            t_2n = pow(t_2, this->order);
        }
        return exp(-xt_2)*t_2n;
    }

protected:
    const int order;
    const double pos;
};

class BoysFunction : public RepresentableFunction<1> {
public:
    BoysFunction(int n, double prec = 1.0e-12)
            : RepresentableFunction<1>(),
              Q(0),
              order(n) {
        int k = 11;
        InterpolatingBasis basis(k);
        BoundingBox<1> world;
        MultiResolutionAnalysis<1> MRA(world, basis);
        this->Q = new MWProjector<1>(MRA, prec);
    }

    virtual ~BoysFunction() {
        if (this->Q != 0) delete this->Q;
    }

    double evalf(const double *r) const {
        int plevel = TelePrompter::getPrintLevel();
        TelePrompter::setPrintLevel(0);
        BoysIntegrand func(this->order, *r);
        FunctionTree<1> *tree = (*this->Q)(func);
        double val = tree->integrate();
        delete tree;
        TelePrompter::setPrintLevel(plevel);
        return val;
    }

protected:
    MWProjector<1> *Q;
    const int order;
};

#endif // BOYSFUNCTION_H

