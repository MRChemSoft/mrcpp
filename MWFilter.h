/*
 *
 *
 *  \date Jul 8, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */

#ifndef MWFILTER_H_
#define MWFILTER_H_

#include <string>
#include <Eigen/Core>

#include "constants.h"
#include "TelePrompter.h"

typedef Eigen::Block<Eigen::MatrixXd, Eigen::Dynamic, Eigen::Dynamic>
FilterBlock;

class MWFilter {
public:
    MWFilter(int k, int t, const std::string &lib = "");
    MWFilter(int t, const Eigen::MatrixXd &data);
    virtual ~MWFilter() { }

    virtual void apply(Eigen::MatrixXd &data) const;
    virtual void applyInverse(Eigen::MatrixXd &data) const;
    virtual void apply(Eigen::VectorXd &data) const;
    virtual void applyInverse(Eigen::VectorXd &data) const;

    int getOrder() const { return this->order; }
    int getType() const { return this->type; }

    const Eigen::MatrixXd &getFilter() const { return this->filter; }
    inline const Eigen::MatrixXd &getSubFilter(int i, int oper = 0) const;
    inline const Eigen::MatrixXd &getCompressionSubFilter(int i) const;
    inline const Eigen::MatrixXd &getReconstructionSubFilter(int i) const;

    static void setDefaultLibrary(const std::string &dir);
    static const std::string &getDefaultLibrary() { return default_filter_lib; }

protected:
    int type;
    int order;
    int dim;

    Eigen::MatrixXd filter; ///< Full MW-transformation matrix
    Eigen::MatrixXd G0;
    Eigen::MatrixXd G1;
    Eigen::MatrixXd H0;
    Eigen::MatrixXd H1;
    // Transpose
    Eigen::MatrixXd G0t;
    Eigen::MatrixXd G1t;
    Eigen::MatrixXd H0t;
    Eigen::MatrixXd H1t;

    std::string H_path;
    std::string G_path;
    static std::string default_filter_lib;

    void setFilterPaths(const std::string &lib);
    void readFilterBin();
    void fillFilterBlocks();
};

const Eigen::MatrixXd& MWFilter::getSubFilter(int i, int oper) const {
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

const Eigen::MatrixXd& MWFilter::getCompressionSubFilter(int i) const {
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

const Eigen::MatrixXd& MWFilter::getReconstructionSubFilter(int i) const {
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

#endif /* MWFILTER_H_ */
