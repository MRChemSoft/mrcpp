/*
 *
 *
 *  \date Jul 8, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */

#include <fstream>

#include "constants.h"
#include "config.h"

#include "MWFilter.h"
#include "utils/Printer.h"
#include "eigen_disable_warnings.h"

using namespace std;
using namespace Eigen;

namespace mrcpp {

string MWFilter::default_filter_lib = MW_FILTER_DIR;

MWFilter::MWFilter(int k, int t, const string &lib)
        : type(t),
          order(k) {
    if (this->order < 1 or this->order > MaxOrder) {
        MSG_FATAL("Invalid filter order: " << this->order);
    }
    switch (this->type) {
    case (Interpol):
    case (Legendre):
        break;
    default:
        MSG_ERROR("Unknown filter type: " << this->type);
    }
    char *ep = getenv("MRCPP_FILTER_DIR");
    if (ep != 0) {
        default_filter_lib = *ep;
    }
    int K = this->order + 1;
    setFilterPaths(lib);

    this->filter = MatrixXd(2*K, 2*K);
    this->filter.setZero();

    readFilterBin();
    fillFilterBlocks();
}

MWFilter::MWFilter(int t, const MatrixXd &data) {
    this->type = t;
    this->order = data.cols()/2 - 1;
    if (this->order < 0 or this->order > MaxOrder) {
        MSG_FATAL("Invalid filter order " << this->order);
    }
    switch (this->type) {
    case (Interpol):
    case (Legendre):
        break;
    default:
        MSG_ERROR("Unknown filter type: " << type);
    }

    this->filter = data;
    fillFilterBlocks();
}

void MWFilter::setDefaultLibrary(const std::string &dir) {
    if (dir.empty()) {
        MSG_ERROR("No directory specified!");
    }
    default_filter_lib = dir;
}

void MWFilter::fillFilterBlocks() {
    int K = this->order + 1;
    this->G0 = this->filter.block(0, 0, K, K);
    this->G1 = this->filter.block(0, K, K, K);
    this->H0 = this->filter.block(K, 0, K, K);
    this->H1 = this->filter.block(K, K, K, K);
    this->G0t = this->G0.transpose();
    this->G1t = this->G1.transpose();
    this->H0t = this->H0.transpose();
    this->H1t = this->H1.transpose();
}

const MatrixXd& MWFilter::getSubFilter(int i, int oper) const {
    switch (oper) {
    case (Compression):
        switch (i) {
        case (0):
            return this->H0t;
        case (1):
            return this->H1t;
        case (2):
            return this->G0t;
        case (3):
            return this->G1t;
        default:
            MSG_FATAL("Filter index out of bounds");
        }
        break;
    case (Reconstruction):
        switch (i) {
        case (0):
            return this->H0;
        case (1):
            return this->G0;
        case (2):
            return this->H1;
        case (3):
            return this->G1;
        default:
            MSG_FATAL("Filter index out of bounds");
        }
        break;
    default:
        MSG_FATAL("Invalid wavelet transformation");
    }
}

const MatrixXd& MWFilter::getCompressionSubFilter(int i) const {
    switch (i) {
    case (0):
        return this->H0t;
    case (1):
        return this->H1t;
    case (2):
        return this->G0t;
    case (3):
        return this->G1t;
    default:
        MSG_FATAL("Filter index out of bounds");
    }
}

const MatrixXd& MWFilter::getReconstructionSubFilter(int i) const {
    switch (i) {
    case (0):
        return this->H0;
    case (1):
        return this->G0;
    case (2):
        return this->H1;
    case (3):
        return this->G1;
    default:
        MSG_FATAL("Filter index out of bounds");
    }
}

void MWFilter::apply(MatrixXd &data) const {
    if (data.rows() != this->filter.cols()) {
        INVALID_ARG_ABORT
    }
    data = this->filter * data;
}

void MWFilter::applyInverse(MatrixXd &data) const {
    if (data.rows() != this->filter.cols()) {
        INVALID_ARG_ABORT
    }
    data = this->filter.transpose() * data;
}

void MWFilter::apply(VectorXd &data) const {
    if (data.rows() != this->filter.cols()) {
        INVALID_ARG_ABORT
    }
    data = this->filter * data;
}

void MWFilter::applyInverse(VectorXd &data) const {
    if (data.rows() != this->filter.cols()) {
        INVALID_ARG_ABORT
    }
    data = this->filter.transpose() * data;
}

void MWFilter::setFilterPaths(const string &lib) {
    ostringstream oss;
    oss << this->order;
    string ordr = oss.str();
    string flib;

    if (lib.empty()) {
        flib = default_filter_lib;
    } else {
        flib = lib;
    }
    switch (this->type) {
    case (Interpol):
        this->H_path = flib + "/I_H0_" + ordr;
        this->G_path = flib + "/I_G0_" + ordr;
        break;
    case (Legendre):
        this->H_path = flib + "/L_H0_" + ordr;
        this->G_path = flib + "/L_G0_" + ordr;
        break;
    default:
        MSG_FATAL("Invalid filter type " << this->type);
    }
}

void MWFilter::readFilterBin() {
    ifstream H_fis(this->H_path.c_str(), ios::binary);
    ifstream G_fis(this->G_path.c_str(), ios::binary);

    if (H_fis.fail()) MSG_FATAL("Could not open filter: " << this->H_path);
    if (G_fis.fail()) MSG_FATAL("Could not open filter: " << this->G_path);

    int K = this->order + 1;
    double dH[K];
    double dG[K];
    int i, j;

    FilterBlock g0 = this->filter.block(0, 0, K, K);
    FilterBlock g1 = this->filter.block(0, K, K, K);
    FilterBlock h0 = this->filter.block(K, 0, K, K);
    FilterBlock h1 = this->filter.block(K, K, K, K);

    /* read H0 and G0 from disk */
    for (i = 0; i < K; i++) {
        H_fis.read((char *) dH, sizeof(double) * K);
        G_fis.read((char *) dG, sizeof(double) * K);
        for (j = 0; j < K; j++) {
            g0(i, j) = dG[j]; // G0
            h0(i, j) = dH[j]; // H0

        }
    }

    /* fill H1 and G1 according to symmetry */
    switch (this->type) {
    case Interpol:
        for (i = 0; i < K; i++) {
            for (j = 0; j < K; j++) {
                g1(i, j) = pow(-1.0, i + K) * g0(i, K - j - 1);
                h1(i, j) = h0(K - i - 1, K - j - 1);
            }
        }
        break;
    case Legendre:
        for (i = 0; i < K; i++) {
            for (j = 0; j < K; j++) {
                g1(i, j) = pow(-1.0, i + j + K) * g0(i, j);
                h1(i, j) = pow(-1.0, i + j) * h0(i, j);

            }
        }
        break;
    }
    G_fis.close();
    H_fis.close();
}

} // namespace mrcpp
