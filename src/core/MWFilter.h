/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2021 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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

#include <string>
#include <vector>

#include <Eigen/Core>

#include "MRCPP/config.h"

namespace mrcpp {

class MWFilter final {
public:
    MWFilter(int k, int t);
    MWFilter(int t, const Eigen::MatrixXd &data);

    void apply(Eigen::MatrixXd &data) const;
    void apply(Eigen::VectorXd &data) const;
    void applyInverse(Eigen::MatrixXd &data) const;
    void applyInverse(Eigen::VectorXd &data) const;

    int getOrder() const { return this->order; }
    int getType() const { return this->type; }

    const Eigen::MatrixXd &getFilter() const { return this->filter; }
    const Eigen::MatrixXd &getSubFilter(int i, int oper = 0) const;
    const Eigen::MatrixXd &getCompressionSubFilter(int i) const;
    const Eigen::MatrixXd &getReconstructionSubFilter(int i) const;

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

private:
    void setFilterPaths(const std::string &lib);
    void fillFilterBlocks();
    void generateBlocks();

    std::string H_path;
    std::string G_path;
};

} // namespace mrcpp
