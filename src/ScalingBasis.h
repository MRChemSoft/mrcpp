#pragma once

#include <vector>

#include "constants.h"

#include "Polynomial.h"

namespace mrcpp {

class ScalingBasis {
public:
    ScalingBasis(int k, int t);
    virtual ~ScalingBasis() { }

    void evalf(const double *r, Eigen::MatrixXd &vals) const;

    Polynomial &getFunc(int k) { return this->funcs[k]; }
    const Polynomial &getFunc(int k) const { return this->funcs[k]; }

    int getScalingType() const { return this->type; }
    int getScalingOrder() const { return this->order; }
    int getQuadratureOrder() const { return this->order + 1; }

    const Eigen::MatrixXd &getQuadratureValues() const { return this->quadVals; }
    const Eigen::MatrixXd &getCVMap(int operation) const;

    bool operator==(const ScalingBasis &basis) const;
    bool operator!=(const ScalingBasis &basis) const;

    friend std::ostream& operator<<(std::ostream &o, const ScalingBasis &bas) {
        o << " polynomial order = " << bas.getScalingOrder() << std::endl;
        if (bas.getScalingType() == Legendre) {
            o << " polynomial type  = Legendre";
        } else if (bas.getScalingType() == Interpol) {
            o << " polynomial type  = Interpolating";
        } else {
            o << " polynomial type  = Unknown";
        }
        return o;
    }

protected:
    const int type;
    const int order;
    Eigen::MatrixXd quadVals;// function values at quadrature pts
    Eigen::MatrixXd cvMap;  // coef-value transformation matrix
    Eigen::MatrixXd vcMap;  // value-coef transformation matrix
    std::vector<Polynomial> funcs;
};

}
