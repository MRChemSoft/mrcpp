#include "ConvolutionOperator.h"
#include "CrossCorrelationGenerator.h"
#include "GridGenerator.h"
#include "MWProjector.h"
#include "OperatorTree.h"
#include "GreensKernel.h"
#include "Gaussian.h"
#include "InterpolatingBasis.h"
#include "LegendreBasis.h"
#include "MathUtils.h"

template<int D>
ConvolutionOperator<D>::ConvolutionOperator(const MultiResolutionAnalysis<D> &mra,
                                            double apply, double build)
        : MWOperator<D>(mra, apply),
          build_prec(build) {
    if (this->build_prec < 0.0) {
        this->build_prec = this->apply_prec;
    }
}

template<int D>
ConvolutionOperator<D>::~ConvolutionOperator() {
    this->clearOperator();
    this->clearKernel();
}

template<int D>
void ConvolutionOperator<D>::initializeOperator(GreensKernel &greens_kernel) {
    MultiResolutionAnalysis<1> *kern_mra = this->MRA.getKernelMRA();
    MultiResolutionAnalysis<2> *oper_mra = this->MRA.getOperatorMRA();

    GridGenerator<1> G(*kern_mra);
    MWProjector<1> Q(*kern_mra, this->build_prec/10.0);
    CrossCorrelationGenerator CC(*oper_mra, this->build_prec);

    for (int i = 0; i < greens_kernel.size(); i++) {
        Gaussian<1> &greens_comp = *greens_kernel[i];
        FunctionTree<1> *kern_comp = new FunctionTree<1>(*kern_mra, MaxAllocNodes);
        G(*kern_comp, greens_comp); //Generate empty grid to hold narrow Gaussian
        Q(*kern_comp, greens_comp); //Project Gaussian starting from the empty grid
        OperatorTree *oper_comp = CC(*kern_comp); //Expand 1D kernel into 2D operator

        this->kernel.push_back(kern_comp);
        this->oper.push_back(oper_comp);
    }
    delete kern_mra;
    delete oper_mra;
}

template<int D>
void ConvolutionOperator<D>::clearKernel() {
    for (int i = 0; i < this->kernel.size(); i++) {
        if (this->kernel[i] != 0) delete this->kernel[i];
    }
    this->kernel.clear();
}

template<int D>
double ConvolutionOperator<D>::calcMinDistance(double epsilon) const {
    int maxScale = this->MRA.getMaxScale();
    return sqrt(epsilon * pow(2.0, -maxScale));
}

template<int D>
double ConvolutionOperator<D>::calcMaxDistance() const {
    const double *lb = this->MRA.getWorldBox().getLowerBounds();
    const double *ub = this->MRA.getWorldBox().getUpperBounds();
    return MathUtils::calcDistance(D, lb, ub);
}

template class ConvolutionOperator<1>;
template class ConvolutionOperator<2>;
template class ConvolutionOperator<3>;
