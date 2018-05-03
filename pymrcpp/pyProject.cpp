/*
 *  \date Apr 26, 2018
 *  \author Magnar Bjørgve <magnar.bjorgve@uit.no> \n
 *          Hylleraas Centre for Quantum Molecular Sciences \n
 *          UiT - The Arctic University of Norway
 */

#include "treebuilders/project.h"
#include "trees/FunctionTree.h"
#include "PyAnalyticFunction.h"

#include "pyProject.h"


void project3D(double prec,
               mrcpp::FunctionTree<3> &out,
               std::function<double (double x, double y, double z)> func,
               int maxIter) {
    PyAnalyticFunction3D inp(func);
    mrcpp::project(prec, out, inp, maxIter);
}
