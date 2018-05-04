/*
 *  \date May 02, 2018
 *  \author Magnar Bjørgve <magnar.bjorgve@uit.no> \n
 *          Hylleraas Centre for Quantum Molecular Sciences \n
 *          UiT - The Arctic University of Norway
 */

#pragma once

#include "functions/RepresentableFunction.h"

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"

class PyAnalyticFunction3D : public mrcpp::RepresentableFunction<3> {
public:
    PyAnalyticFunction3D(std::function<double (double x, double y, double z)> f)
            : func(f) { }

    virtual double evalf(const double *r) const {
        double val = 0.0;
        if (not this->outOfBounds(r)) val = this->func(r[0], r[1], r[2]);
        return val;
    }
protected:
    std::function<double (double x, double y, double z)> func;
};

class PyAnalyticFunction2D : public mrcpp::RepresentableFunction<2> {
public:
    PyAnalyticFunction2D(std::function<double (double x, double y)> f)
            : func(f) { }

    virtual double evalf(const double *r) const {
        double val = 0.0;
        if (not this->outOfBounds(r)) val = this->func(r[0], r[1]);
        return val;
    }
protected:
    std::function<double (double x, double y)> func;
};

class PyAnalyticFunction1D : public mrcpp::RepresentableFunction<1> {
public:
    PyAnalyticFunction1D(std::function<double (double x)> f)
            : func(f) { }

    virtual double evalf(const double *r) const {
        double val = 0.0;
        if (not this->outOfBounds(r)) val = this->func(r[0]);
        return val;
    }
protected:
    std::function<double (double x)> func;
};
