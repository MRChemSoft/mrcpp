#include "ConvolutionOperator.h"
#include "CrossCorrelationGenerator.h"
#include "MWProjector.h"
#include "OperatorTree.h"
#include "GreensKernel.h"
#include "Gaussian.h"
#include "InterpolatingBasis.h"
#include "LegendreBasis.h"

template<int D>
ConvolutionOperator<D>::ConvolutionOperator(const MultiResolutionAnalysis<D> &mra,
                                            double apply, double build, int iter)
        : MWOperator<D>(mra, apply, iter),
          build_prec(build) {
    if (this->build_prec < 0.0) {
        this->build_prec = this->apply_prec;
    }
}

template<int D>
ConvolutionOperator<D>::~ConvolutionOperator() {
    clearKernel();
}

template<int D>
void ConvolutionOperator<D>::initializeOperator(GreensKernel &greens_kernel) {
    double proj_prec = this->build_prec/100.0;
    double ccc_prec = this->build_prec/10.0;

    MultiResolutionAnalysis<1> *kern_mra = this->getKernelMRA();
    MultiResolutionAnalysis<2> *oper_mra = this->getOperatorMRA();

    MWProjector<1> Q(*kern_mra, proj_prec);
    CrossCorrelationGenerator G(*oper_mra, ccc_prec);

    for (int i = 0; i < greens_kernel.size(); i++) {
        Gaussian<1> &greens_comp = *greens_kernel[i];
        FunctionTree<1> *kern_comp = Q(greens_comp);
        OperatorTree *oper_comp = G(*kern_comp);

        this->kernel.push_back(*kern_comp);
        this->oper.push_back(*oper_comp);
    }
    delete kern_mra;
    delete oper_mra;
}

template<int D>
void ConvolutionOperator<D>::clearKernel() {
    for (int i = 0; i < this->kernel.size(); i++) {
        if (this->kernel[i] != 0) delete this->kernel[i];
        this->kernel.clear();
    }
}

template<int D>
MultiResolutionAnalysis<1>* ConvolutionOperator<D>::getKernelMRA() const {
    MultiResolutionAnalysis<1> *mra = 0;
    const BoundingBox<D> &box = this->MRA.getWorldBox();
    const ScalingBasis &basis = this->MRA.getScalingBasis();

    int type = basis.getScalingType();
    int kern_order = 2*basis.getScalingOrder() + 1;

    ScalingBasis *kern_basis = 0;
    if (type == Interpol) {
        kern_basis = new InterpolatingBasis(kern_order);
    } else if (type == Legendre) {
        kern_basis = new LegendreBasis(kern_order);
    } else {
        MSG_ERROR("Invalid scaling type");
        return mra;
    }

    int max_l = 0;
    for (int i = 0; i < D; i++) {
        if (box.size(i) > max_l) {
            max_l = box.size(i);
        }
    }
    int start_l = -max_l;
    int tot_l = 2*max_l;
    NodeIndex<1> idx(box.getScale(), &start_l);
    BoundingBox<1> kern_box(idx, &tot_l);
    mra = new MultiResolutionAnalysis<1>(kern_box, *kern_basis);
    delete kern_basis;
    return mra;
}

template class ConvolutionOperator<1>;
template class ConvolutionOperator<2>;
template class ConvolutionOperator<3>;
