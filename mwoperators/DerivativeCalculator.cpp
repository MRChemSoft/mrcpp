#include "DerivativeCalculator.h"
#include "QuadratureCache.h"
#include "InterpolatingBasis.h"
#include "LegendreBasis.h"
#include "MWNode.h"

using namespace std;
using namespace Eigen;

DerivativeCalculator::DerivativeCalculator(const ScalingBasis &basis,
                                           double a,
                                           double b)
        : A(a), B(b) {
    int kp1 = basis.getQuadratureOrder();
    this->K = MatrixXd::Zero(kp1, kp1);
    this->valueZero = VectorXd::Zero(kp1);
    this->valueOne = VectorXd::Zero(kp1);
    calcKMatrix(basis);
    calcValueVectors(basis);
}

void DerivativeCalculator::calcValueVectors(const ScalingBasis &basis) {
    int kp1 = basis.getQuadratureOrder();
    double sqrtCoef[kp1];
    for (int i = 0; i < kp1; i++) {
        sqrtCoef[i] = sqrt(2.0 * i + 1.0);
    }

    switch (basis.getScalingType()) {
    case Interpol:
        for (int i = 0; i < kp1; i++) {
            this->valueZero(i) = basis.getFunc(i).evalf(0.0);
            this->valueOne(i) = basis.getFunc(i).evalf(1.0);
        }
        break;
    case Legendre:
        for (int i = 0; i < kp1; i++) {
            double val = sqrtCoef[i];
            this->valueOne(i) = val;
            if (IS_ODD(i)) val *= -1.0;
            this->valueZero(i) = val;
        }
        break;
    default:
        MSG_ERROR("Invalid scaling type");
        break;
    }
}

void DerivativeCalculator::calcKMatrix(const ScalingBasis &basis) {
    int kp1 = basis.getQuadratureOrder();
    double sqrtCoef[kp1];
    for (int i = 0; i < kp1; i++) {
        sqrtCoef[i] = sqrt(2.0 * i + 1.0);
    }
    getQuadratureCache(qCache);
    const VectorXd &roots = qCache.getRoots(kp1);
    const VectorXd &weights = qCache.getWeights(kp1);
    VectorXd sqrtWeights = weights.array().sqrt();

    switch (basis.getScalingType()) {
    case Interpol:
        for (int i = 0; i < kp1; i++) {
            Polynomial poly = basis.getFunc(i).calcDerivative();
            for (int j = 0; j < kp1; j++) {
                this->K(i,j) = 2.0 * sqrtWeights(j) * poly.evalf(roots(j));
            }
        }
        break;
    case Legendre:
        for (int i = 0; i < kp1; i++) {
            for (int j = i; j < kp1; j++) {
                if (IS_ODD((i - j))) {
                    this->K(j,i) = 2.0 * sqrtCoef[i] * sqrtCoef[j];
                }
            }
        }
        break;
    default:
        MSG_ERROR("Invalid scaling type");
        break;
    }
}

void DerivativeCalculator::calcNode(MWNode<2> &node) {
    NOT_IMPLEMENTED_ABORT;
    /*
    node.zeroCoefs();
    int kp1 = node.getKp1();
    int kp1_d = node.getKp1_d();
    int l = node.getTranslation()[1] - node.getTranslation()[0];
    double two_np1 = pow(2.0, node.getScale() + 1);

    double a = this->A;
    double b = this->B;
    switch (l) {
    case 1:
        for (int i = 0; i < kp1; i++) {
            double zero_i = this->valueZero(i);
            for (int j = 0; j < kp1; j++) {
                double one_j = this->valueOne(j);
                int idx = i*kp1 + j;
                node.getCoefs()[1*kp1_d + idx] = -b * zero_i * one_j;
            }
        }
        node.getCoefs().segment(1*kp1_d, kp1_d) *= two_np1;
        break;
    case 0:
        for (int i = 0; i < kp1; i++) {
            double zero_i = this->valueZero(i);
            double one_i = this->valueOne(i);
            for (int j = 0; j < kp1; j++) {
                double K_ij = this->K(i,j);
                double one_j = this->valueOne(j);
                double one_ij = one_i * one_j;
                double zero_j = this->valueZero(j);
                double zero_ij = zero_i * zero_j;

                int idx = i*kp1 + j;
                node.getCoefs()[0*kp1_d + idx] = (1.0-a)*one_ij - (1.0-b)*zero_ij - K_ij;
                node.getCoefs()[1*kp1_d + idx] =  a * one_i * zero_j;
                node.getCoefs()[2*kp1_d + idx] = -b * zero_i * one_j;
                node.getCoefs()[3*kp1_d + idx] = (1.0-a)*one_ij - (1.0-b)*zero_ij - K_ij;
            }
        }
        node.getCoefs() *= two_np1;
        break;
    case -1:
        for (int i = 0; i < kp1; i++) {
            double one_i = this->valueOne(i);
            for (int j = 0; j < kp1; j++) {
                double zero_j = this->valueZero(j);
                int idx = i*kp1 + j;
                node.getCoefs()[2*kp1_d + idx] = a * one_i * zero_j;
            }
        }
        node.getCoefs().segment(2*kp1_d, kp1_d) *= two_np1;
        break;
    default:
        MSG_ERROR("This translation should not occour");
        break;
    }
    node.mwTransform(Compression);
    node.setHasCoefs();
    node.calcNorms();
    */
}
