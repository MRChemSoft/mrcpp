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

/**
 *  Basic class for representing functions in a multiwavelet
 *  representation.
 */

#pragma once

#include "MWTree.h"

namespace mrcpp {

template <int D> class FunctionTree final : public MWTree<D> {
public:
    FunctionTree(const MultiResolutionAnalysis<D> &mra, SharedMemory *sh_mem = nullptr);
    FunctionTree(const FunctionTree<D> &tree) = delete;
    FunctionTree<D> &operator=(const FunctionTree<D> &tree) = delete;
    ~FunctionTree();

    void clear();

    double integrate() const;
    double evalf(const Coord<D> &r);

    void getEndValues(Eigen::VectorXd &data);
    void setEndValues(Eigen::VectorXd &data);

    void saveTree(const std::string &file);
    void loadTree(const std::string &file);

    // In place operations
    void square();
    void power(double p);
    void rescale(double c);
    void normalize();
    void add(double c, FunctionTree<D> &inp);
    void multiply(double c, FunctionTree<D> &inp);

    int getNChunks();
    int getNChunksUsed();

    int crop(double prec, double splitFac = 1.0, bool absPrec = true);

    FunctionNode<D> &getEndFuncNode(int i) { return static_cast<FunctionNode<D> &>(this->getEndMWNode(i)); }
    FunctionNode<D> &getRootFuncNode(int i) { return static_cast<FunctionNode<D> &>(this->rootBox.getNode(i)); }

    SerialFunctionTree<D> *getSerialFunctionTree() { return static_cast<SerialFunctionTree<D> *>(this->serialTree_p); }
    void printSerialIndices();

    const FunctionNode<D> &getEndFuncNode(int i) const {
        return static_cast<const FunctionNode<D> &>(this->getEndMWNode(i));
    }
    const FunctionNode<D> &getRootFuncNode(int i) const {
        return static_cast<const FunctionNode<D> &>(this->rootBox.getNode(i));
    }

protected:
    std::ostream &print(std::ostream &o);
};

} // namespace mrcpp
