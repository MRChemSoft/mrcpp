/*
 *  \date Apr 26, 2018
 *  \author Magnar Bjørgve <magnar.bjorgve@uit.no> \n
 *          Hylleraas Centre for Quantum Molecular Sciences \n
 *          UiT - The Arctic University of Norway
 */

#pragma once

#include <functional>

#include "mrcpp_declarations.h"


void project3D(double prec, mrcpp::FunctionTree<3> &out,
               std::function<double (double x, double y, double z)> func, int maxIter = -1);

void project2D(double prec, mrcpp::FunctionTree<2> &out,
               std::function<double (double x, double y)> func, int maxIter = -1);

void project1D(double prec, mrcpp::FunctionTree<1> &out,
               std::function<double (double x)> func, int maxIter = -1);
