/*
 * 
 *
 *  \date Jul 7, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */

#ifndef HELMHOLTZKERNEL_H_
#define HELMHOLTZKERNEL_H_

#include "GreensKernel.h"

class HelmholtzKernel: public GreensKernel {
public:
	template<int D>
	HelmholtzKernel(double _mu, const MultiResolutionAnalysis<D> &mra,
					bool _delta = false, double _deltaMu = 1.0):
			GreensKernel(mra) {
		mu = _mu;
		delta = _delta;
		deltaMu = _deltaMu;
		compKernRepr();
	}
	virtual ~HelmholtzKernel() {}
	// TODO Auto-generated destructor stub

	double getMu() const {
		return mu;
	}
protected:
	double mu; /**< exponent */
	bool delta; /**< Delta or full operator */
	bool deltaMu; /**< Delta or full operator */
	virtual void compKernRepr();
};

#endif /* HELMHOLTZKERNEL_H_ */
