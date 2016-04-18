/*
 * \breif
 */

#ifndef GREENSKERNEL_H_
#define GREENSKERNEL_H_

#include "Gaussian.h"
#include "GaussExp.h"
#include "GaussFunc.h"
#include "constants.h"
#include "MultiResolutionAnalysis.h"

const int MAX_SEP_RANK = 200;
const double MIN_DIST = 1.0e-14;
const double MAX_DIST = 999.0;
const double MIN_EPSILON = 1.e-14;
const double MAX_EPSILON = 1.e-1;

class GreensKernel {
public:
    template <int D> GreensKernel(const MultiResolutionAnalysis<D> &mra) {
        setDefaultParams();
        genKernBox(mra.getWorldBox());
    }

    virtual ~GreensKernel() { }

    virtual double evalf(const double *r) { return kern.evalf(r); }

    GaussExp<1> apply(Gaussian<1> &f);
    GaussExp<1> apply(GaussExp<1> &g);
    GaussFunc<1> apply1d(Gaussian<1> &f);
    GaussExp<1> apply1d(GaussExp<1> &g);

    Gaussian<1> &getKernComp(int n) { return kern.getFunc(n); }
    GaussExp<1> &getKern() { return kern; }
    const BoundingBox<1> &getKernBox() const { return kernBox; }

    int getNTerms() const { return kern.size(); }
    int getCompMin() const { return compMin; }
    int getCompMax() const { return compMax; }
    double getRMin() const { return rMin; }
    double getRMax() const { return rMax; }
    double getExpPrec() const { return expPrec; }
    double getProjPrec() const { return projPrec; }

    void setRMin(double r, bool reset = true);
    void setRMax(double r, bool reset = true);
    void setCompMin(int n, bool reset = true);
    void setCompMax(int n, bool reset = true);
    void setExpPrec(double prec, bool reset = true);
    void setProjPrec(double prec, bool reset = true);
    void setKernelParam(double expPrec = -1.0, double projPrec = -1.0,
    double rMin = -1.0,	double rMax = -1.0,	int nMin = -1, int nMax = -1);

    friend std::ostream& operator <<(std::ostream &o, const GreensKernel &kernel) {
        o << "Kernel: " << std::endl;
        o << "compMin:  " << kernel.compMin << std::endl;
        o << "compMax:  " << kernel.compMax << std::endl;
        o << "rMin:     " << kernel.rMin << std::endl;
        o << "rMax:     " << kernel.rMax << std::endl;
        o << "expPrec:  " << kernel.expPrec << std::endl;
        o << "projPrec: " << kernel.projPrec << std::endl << std::endl;
        o << kernel.kern << std::endl;
        return o;
    }
    static void setDefaultRMin(double rMin);
    static void setDefaultRMax(double rMax);
    static void setDefaultCompMin(int compMin);
    static void setDefaultCompMax(int compMax);
    static void setDefaultExpPrec(double expPrec);
    static void setDefaultProjPrec(double projPrec);
protected:
    static double defaultRMin;
    static double defaultRMax;
    static int defaultCompMin;
    static int defaultCompMax;
    static double defaultExpPrec;
    static double defaultProjPrec;

    int compMin; /**< lowest comp used */
    int compMax; /**< highest comp used */
    double rMin; /**< lower extreme */
    double rMax; /**< upper extreme  (not used) */
    double expPrec; /**< precision of kernel expansion*/
    double projPrec; /**< precision of PROJECTION of kernel expansion*/
    BoundingBox<1> kernBox;

    virtual void compKernRepr() = 0;
    GaussExp<1> kern; //< Kernel representation as GaussExp
    GaussExp<1> convolute(const Gaussian<1> &g);
    GaussFunc<1> convolute1(int k, const Gaussian<1> &g);

    void setDefaultParams() {
        compMin = defaultCompMin;
        compMax = defaultCompMax;
        rMin = defaultRMin;
        rMax = defaultRMax;
        expPrec = defaultExpPrec;
        projPrec = defaultProjPrec;
    }
    template<int D>
    void genKernBox(const BoundingBox<D> &box) {
        int transl = 0;
        int maxn = 0;
        for (int i = 0; i < D; i++) {
            if (box.size(i) > maxn) {
                maxn = box.size(i);
                transl = box.getCornerIndex().getTranslation(i);
            }
        }
        int l = transl - maxn;
        maxn *= 2;
        NodeIndex<1> idx(box.getScale(), &l);
        kernBox = BoundingBox<1>(idx, &maxn);
    }
};

#endif /* GREENSKERNEL_H_ */
