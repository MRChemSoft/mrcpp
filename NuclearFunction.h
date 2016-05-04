#ifndef NUCLEARFUNCTION_H
#define NUCLEARFUNCTION_H

#include <vector>

#include "RepresentableFunction.h"
#include "MathUtils.h"

class NuclearFunction : public RepresentableFunction<3> {
public:
    NuclearFunction() : RepresentableFunction<3>() {
        this->constant = -1.0/(3.0*sqrt(pi));
    }

    virtual ~NuclearFunction() { }

    void addNucleus(const double pos[3], double Z, double smooth) {
        this->x_coords.push_back(pos[0]);
        this->y_coords.push_back(pos[1]);
        this->z_coords.push_back(pos[2]);
        this->charges.push_back(Z);
        this->smoothParam.push_back(smooth);
    }

    double evalf(const double *x) const {
        double result = 0.0;
        for (int i = 0; i < this->x_coords.size(); i++) {
            double xyz[3];
            xyz[0] = this->x_coords[i];
            xyz[1] = this->y_coords[i];
            xyz[2] = this->z_coords[i];
            double r1 = MathUtils::calcDistance(3, xyz, x);
            r1 *= 1.0/this->smoothParam[i];
            double r2 = pow(r1, 2.0);
            double partResult = -erf(r1)/r1;
            partResult += this->constant*(exp(-r2) + 16.0*exp(-4.0*r2));
            result += this->charges[i]*partResult/this->smoothParam[i];
        }
        return result;
    }

    bool isVisibleAtScale(int scale, int nQuadPts) const {
        double minSmooth = this->smoothParam[0];
        for (int i = 1; i < this->smoothParam.size(); i++) {
            if (this->smoothParam[i] < minSmooth) {
                minSmooth = this->smoothParam[i];
            }
        }
        double stdDeviation = pow(minSmooth, -0.5);
        int visibleScale = int (floor(log2(nQuadPts*5.0*stdDeviation)));
        if (scale < visibleScale) {
            return false;
        }
        return true;
    }

    bool isZeroOnInterval(const double *a, const double *b) const {
        int nNucs = this->x_coords.size();
        int totSplit = 0;
        for (int i = 0; i < nNucs; i++) {
            double x_i = this->x_coords[i];
            double y_i = this->y_coords[i];
            double z_i = this->z_coords[i];
            int split = 1;
            if (a[0] > x_i or b[0] < x_i) split = 0;
            if (a[1] > y_i or b[1] < y_i) split = 0;
            if (a[2] > z_i or b[2] < z_i) split = 0;
            totSplit += split;
        }
        if (totSplit == 0) {
            return true;
        } else {
            return false;
        }
    }

protected:
    double constant;
    std::vector<double> smoothParam;
    std::vector<double> charges;
    std::vector<double> x_coords;
    std::vector<double> y_coords;
    std::vector<double> z_coords;
};

#endif // NUCLEARFUNCTION_H
