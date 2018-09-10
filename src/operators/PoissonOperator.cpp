#include "PoissonOperator.h"
#include "PoissonKernel.h"
#include "utils/Printer.h"

namespace mrcpp {

PoissonOperator::PoissonOperator(const MultiResolutionAnalysis<3> &mra, double pr)
        : ConvolutionOperator<3>(mra, pr) {
    int oldlevel = Printer::setPrintLevel(0);
    double epsilon = this->prec/10.0;
    double r_min = calcMinDistance(mra, epsilon);
    double r_max = calcMaxDistance(mra);
    PoissonKernel poisson_kernel(epsilon, r_min, r_max);
    // Rescale for application in 3D
    poisson_kernel.rescale(3);
    initializeOperator(poisson_kernel);
    Printer::setPrintLevel(oldlevel);
}

} // namespace mrcpp
