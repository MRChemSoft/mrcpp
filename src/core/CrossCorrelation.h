/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2019 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
 *
 * This file is part of MRCPP.
 *
 * MRCPP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRCPP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRCPP.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRCPP, see:
 * <https://mrcpp.readthedocs.io/>
 */

#pragma once

#pragma GCC system_header
#include <Eigen/Core>

#include <string>

namespace mrcpp {

class CrossCorrelation final {
public:
    CrossCorrelation(int k, int t, const std::string &lib = "");
    CrossCorrelation(int t, const Eigen::MatrixXd &ldata, const Eigen::MatrixXd &rdata);

    int getType() const { return this->type; }
    int getOrder() const { return this->order; }
    const Eigen::MatrixXd &getLMatrix() const { return this->Left; }
    const Eigen::MatrixXd &getRMatrix() const { return this->Right; }

    static void setDefaultLibrary(const std::string &dir);
    static const std::string &getDefaultLibrary() { return default_ccc_lib; }

protected:
    int type;
    int order;

    Eigen::MatrixXd Left;
    Eigen::MatrixXd Right;

    std::string L_path;
    std::string R_path;
    static std::string default_ccc_lib;

    void setCCCPaths(const std::string &lib);
    void readCCCBin();
};

} // namespace mrcpp
