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

#include <map>
#include <memory>

#include "MWTree.h"
#include "NodeAllocator.h"

namespace mrcpp {

/** @class FunctionTree
 *
 * @brief Function representation in MW basis
 *
 * @details
 * Constructing a full grown FunctionTree involves a number of steps,
 * including setting up a memory allocator, constructing root nodes according
 * to the given MRA, building an adaptive tree structure and computing MW
 * coefficients. The FunctionTree constructor does only half of these steps:
 * It takes an MRA argument, which defines the computational domain and scaling
 * basis (these are fixed parameters that cannot be changed after construction).
 * The tree is initialized with a memory allocator and a set of root nodes, but
 * it does not compute any coefficients and the function is initially
 * *undefined*. An undefined FunctionTree will have a well defined tree
 * structure (at the very least the root nodes of the given MRA, but possibly
 * with additional refinement) and its MW coefficient will be allocated but
 * uninitialized, and its square norm will be negative (minus one).
 */

template <int D> class FunctionTree final : public MWTree<D>, public RepresentableFunction<D> {
public:
    FunctionTree(const MultiResolutionAnalysis<D> &mra, const std::string &name)
            : FunctionTree(mra, nullptr, name) {}
    FunctionTree(const MultiResolutionAnalysis<D> &mra, SharedMemory *sh_mem = nullptr, const std::string &name = "nn");
    FunctionTree(const FunctionTree<D> &tree) = delete;
    FunctionTree<D> &operator=(const FunctionTree<D> &tree) = delete;
    ~FunctionTree() override;

    double integrate() const;
    double integrateEndNodes(RepresentableFunction_M &f);
    double evalf_precise(const Coord<D> &r);
    double evalf(const Coord<D> &r) const override;

    int getNGenNodes() const { return getGenNodeAllocator().getNNodes(); }

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
    void absadd(double c, FunctionTree<D> &inp);
    void multiply(double c, FunctionTree<D> &inp);
    void map(FMap fmap);

    int getNChunks() { return this->getNodeAllocator().getNChunks(); }
    int getNChunksUsed() { return this->getNodeAllocator().getNChunksUsed(); }

    int crop(double prec, double splitFac = 1.0, bool absPrec = true);

    FunctionNode<D> &getEndFuncNode(int i) { return static_cast<FunctionNode<D> &>(this->getEndMWNode(i)); }
    FunctionNode<D> &getRootFuncNode(int i) { return static_cast<FunctionNode<D> &>(this->rootBox.getNode(i)); }

    NodeAllocator<D> &getGenNodeAllocator() { return *this->genNodeAllocator_p; }
    const NodeAllocator<D> &getGenNodeAllocator() const { return *this->genNodeAllocator_p; }

    const FunctionNode<D> &getEndFuncNode(int i) const { return static_cast<const FunctionNode<D> &>(this->getEndMWNode(i)); }
    const FunctionNode<D> &getRootFuncNode(int i) const { return static_cast<const FunctionNode<D> &>(this->rootBox.getNode(i)); }

    void deleteGenerated();
    void deleteGeneratedParents();

    void makeCoeffVector(std::vector<double *> &coefs,
                         std::vector<int> &indices,
                         std::vector<int> &parent_indices,
                         std::vector<double> &scalefac,
                         int &max_index,
                         MWTree<D> &refTree,
                         std::vector<MWNode<D> *> *refNodes = nullptr);
    void makeTreefromCoeff(MWTree<D> &refTree, std::vector<double *> coefpVec, std::map<int, int> &ix2coef, double absPrec, const std::string &mode = "adaptive");
    void appendTreeNoCoeff(MWTree<D> &inTree);

    // tools for use of local (nodes are stored in Bank) representation
    int saveNodesAndRmCoeff(); // put all nodes coefficients in Bank and delete all coefficients
protected:
    std::unique_ptr<NodeAllocator<D>> genNodeAllocator_p{nullptr};
    std::ostream &print(std::ostream &o) const override;

    void allocRootNodes();
};

} // namespace mrcpp
