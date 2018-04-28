#ifndef FACTORY_FUNCTIONS_H
#define FACTORY_FUNCTIONS_H

#include "trees/BoundingBox.h"
#include "trees/NodeIndex.h"
#include "trees/MultiResolutionAnalysis.h"
#include "core/InterpolatingBasis.h"
#include "core/LegendreBasis.h"
#include "trees/FunctionTree.h"
#include "functions/GaussFunc.h"
#include "mwutils/Printer.h"

template<class T> void finalize(T **obj) {
    if (obj == 0) MSG_FATAL("Invalid argument");
    if (*obj == 0) MSG_FATAL("Invalid argument");
    delete *obj;
    *obj = 0;
}

template<int D> void initialize(mrcpp::NodeIndex<D> **idx) {
    if (idx == 0) MSG_FATAL("Invalid argument");
    if (*idx != 0) MSG_FATAL("Invalid argument");
    int scale = 1;
    int l[D];
    for (int d = 0; d < D; d++) {
        l[d] = d-1;
    }
    *idx = new mrcpp::NodeIndex<D>(scale, l);
}

template<int D> void testInitial(const mrcpp::NodeIndex<D> *idx) {
    if (idx == 0) MSG_FATAL("Invalid argument");

    const int scale = 1;
    REQUIRE( (scale == idx->getScale()) );

    for (int d = 0; d < D; d++) {
        const int l = d-1;
        REQUIRE( (l == idx->getTranslation(d)) );
        REQUIRE( (l == idx->getTranslation()[d]) );
    }
}

template<int D> void initialize(mrcpp::BoundingBox<D> **box) {
    if (box == 0) MSG_FATAL("Invalid argument");
    if (*box != 0) MSG_FATAL("Invalid argument");

    int nb[D];
    for (int d = 0; d < D; d++) {
        nb[d] = d + 1;
    }
    mrcpp::NodeIndex<D> *nIdx = 0;
    initialize(&nIdx);

    *box = new mrcpp::BoundingBox<D>(*nIdx, nb);
    finalize(&nIdx);
}

template<int D> void testInitial(const mrcpp::BoundingBox<D> *box) {
    if (box == 0) MSG_FATAL("Invalid argument");

    const mrcpp::NodeIndex<D> &cIdx = box->getCornerIndex();
    testInitial<D>(&cIdx);

    REQUIRE( (box->getUnitLength() > 0.0) );

    for (int d = 0; d < D; d++) {
        REQUIRE( (box->getBoxLength(d) > 0.0) );
        REQUIRE( (box->getBoxLengths()[d] > 0.0) );

        REQUIRE( (box->getLowerBound(d) < box->getUpperBound(d)) );
        REQUIRE( (box->getLowerBounds()[d] < box->getUpperBounds()[d]) );
    }

    int tot_boxes = 1;
    for (int d = 0; d < D; d++) {
        const int nb = d + 1;
        REQUIRE( (box->size(d) == nb) );
        tot_boxes *= box->size(d);
    }
    REQUIRE( (box->size() == tot_boxes) );
}

template<int D> void initialize(mrcpp::MultiResolutionAnalysis<D> **mra) {
    if (mra == 0) MSG_FATAL("Invalid argument");
    if (*mra != 0) MSG_FATAL("Invalid argument");

    int k = 5;
    mrcpp::InterpolatingBasis basis(k);
    mrcpp::BoundingBox<D> *world = 0;
    initialize(&world);
    *mra = new mrcpp::MultiResolutionAnalysis<D>(*world, basis);
    finalize(&world);
}

/* Initializing a D-dimensional Gaussian of unit charge */
template<int D> void initialize(mrcpp::GaussFunc<D> **func) {
    double beta = 1.0e4;
    double alpha = std::pow(beta/mrcpp::pi, D/2.0);
    double pos[3] = {-0.2, 0.5, 1.0};
    int pow[3] = {0, 0, 0};
    *func = new mrcpp::GaussFunc<D>(beta, alpha, pos, pow);
}

#endif //FACTORY_FUNCTIONS_H
