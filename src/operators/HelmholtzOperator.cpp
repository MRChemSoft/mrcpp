#include "HelmholtzOperator.h"
#include "HelmholtzKernel.h"
#include "utils/Printer.h"

namespace mrcpp {

HelmholtzOperator::HelmholtzOperator(const MultiResolutionAnalysis<3> &mra,
                                     double m, double pr)
        : ConvolutionOperator<3>(mra, pr), mu(m) {
    int oldlevel = Printer::setPrintLevel(0);
    double epsilon = this->prec/10.0;
    double r_min = calcMinDistance(mra, epsilon);
    double r_max = calcMaxDistance(mra);
    HelmholtzKernel helmholtz_kernel(this->mu, epsilon, r_min, r_max);
    // Rescale for application in 3D
    helmholtz_kernel.rescale(3);
    initializeOperator(helmholtz_kernel);
    Printer::setPrintLevel(oldlevel);
}

} // namespace mrcpp
