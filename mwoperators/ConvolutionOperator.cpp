#include "ConvolutionOperator.h"
#include "CrossCorrelationGenerator.h"
#include "GridGenerator.h"
#include "MWProjector.h"
#include "OperatorTree.h"
#include "GreensKernel.h"
#include "Gaussian.h"
#include "MathUtils.h"

template<int D>
ConvolutionOperator<D>::~ConvolutionOperator() {
    this->clearOperator();
    this->clearKernel();
}

template<int D>
void ConvolutionOperator<D>::initializeOperator(GreensKernel &greens_kernel) {
    GridGenerator<1> G;
    MWProjector<1> Q(this->prec/10.0);
    CrossCorrelationGenerator CC(this->prec);

    for (int i = 0; i < greens_kernel.size(); i++) {
        Gaussian<1> &greens_comp = *greens_kernel[i];
        FunctionTree<1> *kern_comp = new FunctionTree<1>(this->kern_mra, MaxAllocNodes1D);

        G(*kern_comp, greens_comp); //Generate empty grid to hold narrow Gaussian
        Q(*kern_comp, greens_comp); //Project Gaussian starting from the empty grid

        OperatorTree *oper_comp = new OperatorTree(this->oper_mra, this->prec);
        CC(*oper_comp, *kern_comp); //Expand 1D kernel into 2D operator

        this->kernel_exp.push_back(kern_comp);
        this->oper_exp.push_back(oper_comp);
    }
}

template<int D>
void ConvolutionOperator<D>::clearKernel() {
    for (int i = 0; i < this->kernel_exp.size(); i++) {
        if (this->kernel_exp[i] != 0) delete this->kernel_exp[i];
    }
    this->kernel_exp.clear();
}

template<int D>
double ConvolutionOperator<D>::calcMinDistance(const MultiResolutionAnalysis<D> &MRA, double epsilon) const {
    int maxScale = MRA.getMaxScale();
    return sqrt(epsilon * pow(2.0, -maxScale));
}

template<int D>
double ConvolutionOperator<D>::calcMaxDistance(const MultiResolutionAnalysis<D> &MRA) const {
    const double *lb = MRA.getWorldBox().getLowerBounds();
    const double *ub = MRA.getWorldBox().getUpperBounds();
    return MathUtils::calcDistance(D, lb, ub);
}

template class ConvolutionOperator<1>;
template class ConvolutionOperator<2>;
template class ConvolutionOperator<3>;
