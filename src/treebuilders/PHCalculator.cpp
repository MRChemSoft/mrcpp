#include <fstream>

#include "config.h"

#include "PHCalculator.h"
#include "utils/Printer.h"

using Eigen::MatrixXd;

namespace mrcpp {

PHCalculator::PHCalculator(const ScalingBasis &basis, int n)
        : diff_order(n) {
    if (this->diff_order <= 0) NOT_IMPLEMENTED_ABORT;
    if (this->diff_order == 1) readSMatrix(basis, '1');
    if (this->diff_order == 2) readSMatrix(basis, '2');
    if (this->diff_order >= 3) NOT_IMPLEMENTED_ABORT;
}

void PHCalculator::readSMatrix(const ScalingBasis &basis, char n) {
    std::string file;
    std::string path = MW_FILTER_DIR;
    if (basis.getScalingType() == Legendre) file = path + "/L_ph_deriv_" + n + ".txt";
    if (basis.getScalingType() == Interpol) file = path + "/I_ph_deriv_" + n + ".txt";
    if (basis.getScalingOrder() <  0) MSG_FATAL("Scaling order not supported");
    if (basis.getScalingOrder() > 29) MSG_FATAL("Scaling order not supported");

    std::fstream ifs;
    ifs.open(file.c_str());
    if (not ifs) MSG_ERROR("Failed to open file: " << file);
    for (int kp1 = 2; kp1 < 30; kp1++) {
        std::string line;
        getline(ifs, line);
        std::istringstream iss(line);

        int order;
        iss >> order;
        if (order != kp1) MSG_FATAL("Orders no not match");

        MatrixXd data = MatrixXd::Zero(3*kp1, kp1);
        for (int i = 0; i < 3*kp1; i++) {
            getline(ifs, line);
            std::istringstream iss(line);
            for (int j = 0; j < kp1; j++) {
                iss >> data(i,j);
            }
        }
        if (kp1 == (basis.getScalingOrder() + 1)) {
            this->S_p1 = data.block(0*kp1, 0, kp1, kp1);
            this->S_0  = data.block(1*kp1, 0, kp1, kp1);
            this->S_m1 = data.block(2*kp1, 0, kp1, kp1);
            break;
        }
    }
}

void PHCalculator::calcNode(MWNode<2> &node) {
    node.zeroCoefs();
    int np1 = node.getScale() + 1;
    int kp1 = node.getKp1();
    int kp1_d = node.getKp1_d();
    int l = node.getTranslation()[1] - node.getTranslation()[0];
    double two_np1 = std::pow(2.0, this->diff_order*np1);
    double *coefs = node.getCoefs();

    switch (l) {
    case 1:
        for (int i = 0; i < kp1; i++) {
            for (int j = 0; j < kp1; j++) {
                int idx = i*kp1 + j;
                coefs[1*kp1_d+idx] = two_np1*this->S_p1(i,j);
            }
        }
        break;
    case 0:
        for (int i = 0; i < kp1; i++) {
            for (int j = 0; j < kp1; j++) {
                int idx = i*kp1 + j;
                coefs[0*kp1_d+idx] = two_np1*this->S_0(i,j);
                coefs[1*kp1_d+idx] = two_np1*this->S_m1(i,j);
                coefs[2*kp1_d+idx] = two_np1*this->S_p1(i,j);
                coefs[3*kp1_d+idx] = two_np1*this->S_0(i,j);
            }
        }
        break;
    case -1:
        for (int i = 0; i < kp1; i++) {
            for (int j = 0; j < kp1; j++) {
                int idx = i*kp1 + j;
                coefs[2*kp1_d+idx] = two_np1*this->S_m1(i,j);
            }
        }
        break;
    default:
        MSG_ERROR("This translation should not occour");
        break;
    }
    node.mwTransform(Compression);
    node.setHasCoefs();
    node.calcNorms();
}

} // namespace mrcpp
