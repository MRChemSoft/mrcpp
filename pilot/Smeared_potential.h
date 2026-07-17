#pragma once



#ifndef SMEARED_POTENTIAL_H
#define SMEARED_POTENTIAL_H

#include <cmath>
#include <array>
#include <functional>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <numeric>




#include "../api/Printer"
#include "../api/Timer"
#include "../api/Gaussians"
#include "../api/Plotter"
#include "../api/MWFunctions"
#include "../api/MWOperators"




double smear_pot(double radius, int charge, double precision) {

    double r = radius;
    double C = 0.00435 * precision / (charge * charge * charge * charge * charge);
    double factor = std::pow(C, 1.0 / 3.0);
    r = r / factor;


    double first_term = std::erf(r);
    first_term = first_term / r;
    double r_srqd = r * r;

    double second_term_prefactor = (1.0 / (3.0 * std::sqrt(M_PI)));
    double second_term = std::exp(-r_srqd) + 16.0 * std::exp(-4.0 * r_srqd);
    second_term = second_term_prefactor * second_term;

    return (first_term + second_term) * charge / factor;
}



class SmearedPotential {
public:
    static double coulomb_HFYGB(const double& radius, 
                                double charge, 
                                double precision) {
        // Compute squared distance
        double distance = radius;

        // Smoothing function
        auto smoothing_HFYGB = [](double charge, double prec) {
            double factor = 0.00435 * prec / std::pow(charge, 5);
            return std::pow(factor, 1.0 / 3.0);
        };

        // uHFYGB function
        auto uHFYGB = [](double r) {
            double erf_term = std::erf(r) / r;
            double exp_term = (1.0 / (3.0 * std::sqrt(M_PI))) * 
                              (std::exp(-r * r) + 16.0 * std::exp(-4.0 * r * r));
            return erf_term + exp_term;
        };

        double factor = smoothing_HFYGB(charge, precision);
        double value = uHFYGB(distance / factor);

        return charge * value / factor;
    }
};

#endif // SMEARED_POTENTIAL_H