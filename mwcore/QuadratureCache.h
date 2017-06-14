/*
 *
 *
 *  \date Jul 26, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */

#ifndef QUADRATURECACHE_H_
#define QUADRATURECACHE_H_

#pragma GCC system_header
#include <Eigen/Core>

#include "GaussQuadrature.h"
#include "ObjectCache.h"

#define getQuadratureCache(X) QuadratureCache &X=QuadratureCache::getInstance()

class QuadratureCache: public ObjectCache<GaussQuadrature> {
public:
    static QuadratureCache &getInstance() {
        static QuadratureCache theQuadratureCache;
        return theQuadratureCache;
    }

    void load(int order);
    GaussQuadrature &get(int order);

    const Eigen::VectorXd &getRoots(int i) { return get(i).getRoots(); }
    const Eigen::VectorXd &getWeights(int i) { return get(i).getWeights(); }

    void setIntervals(int i);
    void setBounds(double a, double b);

    int getIntervals() const { return this->intervals; }
    double getUpperBound() const { return this->B; }
    double getLowerBound() const { return this->A; }
private:
    double A;
    double B;
    int intervals;

    QuadratureCache();
    virtual ~QuadratureCache();

    QuadratureCache(QuadratureCache const &qc) : ObjectCache<GaussQuadrature>(qc) { }
    QuadratureCache &operator=(QuadratureCache const&) { return *this; }
};

#endif /* QUADRATURECACHE_H_ */
