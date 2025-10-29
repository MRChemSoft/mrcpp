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

/**
 * @file FunctionTree.h
 * @brief Declaration of the FunctionTree class template.
 *
 * @details
 * A FunctionTree represents a scalar field on a multiresolution (MW) grid.
 * It owns the MW-node topology, coefficient storage, and basic utilities
 * for evaluation, integration, and in–place algebra on the represented
 * function. Construction initializes the tree structure (root nodes and
 * allocator) according to a given MultiResolutionAnalysis (MRA), but does
 * not compute coefficients; initially the function is undefined and the
 * tree's square norm is negative to signal this state.
 */

/**
 * @class FunctionTree
 * @tparam D Spatial dimension (1, 2, or 3).
 * @tparam T Scalar type (double or ComplexDouble).
 * @brief Function representation in the MW basis with adaptive topology.
 *
 * @details
 * The class derives from MWTree (topology and node management) and
 * RepresentableFunction (evaluation interface). Typical workflows build
 * or refine the tree via calculators/adaptors and then apply algebraic
 * transforms in place.
 */
template <int D, typename T>
class FunctionTree final : public MWTree<D, T>, public RepresentableFunction<D, T> {
public:
    /**
     * @brief Construct a tree bound to an MRA with a user label.
     * @param mra Multi-resolution analysis (domain and basis).
     * @param name Optional textual name of the function.
     *
     * @note Coefficients are not computed by the constructor.
     */
    FunctionTree(const MultiResolutionAnalysis<D> &mra, const std::string &name)
            : FunctionTree(mra, nullptr, name) {}

    /**
     * @brief Construct a tree bound to an MRA with optional shared memory and name.
     * @param mra Multi-resolution analysis (domain and basis).
     * @param sh_mem Optional shared-memory arena for coefficient storage.
     * @param name Optional textual name of the function.
     *
     * @details Root nodes and allocators are created. The function is
     * undefined until coefficients are computed by a builder/calculator.
     */
    FunctionTree(const MultiResolutionAnalysis<D> &mra,
                 SharedMemory<T> *sh_mem = nullptr,
                 const std::string &name = "nn");

    /// Deleted copy semantics (trees are heavy objects).
    FunctionTree(const FunctionTree<D, T> &tree) = delete;
    FunctionTree<D, T> &operator=(const FunctionTree<D, T> &tree) = delete;

    /// Virtual destructor.
    ~FunctionTree() override;

    /**
     * @brief Integrate the represented function over the world domain.
     * @return Integral value.
     */
    T integrate() const;

    /**
     * @brief Integrate only end nodes against a provided analytic function.
     * @param f RepresentableFunction used as integrand partner.
     * @return Integral value as double.
     *
     * @details Useful for quadrature-like post-processing on the current grid.
     */
    double integrateEndNodes(RepresentableFunction_M &f);

    /**
     * @brief Evaluate with high accuracy at a given coordinate.
     * @param r Physical coordinate.
     * @return Function value.
     *
     * @details May be more expensive than evalf due to stricter handling.
     */
    T evalf_precise(const Coord<D> &r);

    /**
     * @brief Evaluate the function at a given coordinate.
     * @param r Physical coordinate.
     * @return Function value.
     */
    T evalf(const Coord<D> &r) const override;

    /**
     * @brief Number of generated (non-root) nodes currently alive.
     * @return Count of nodes managed by the generated-node allocator.
     */
    int getNGenNodes() const { return getGenNodeAllocator().getNNodes(); }

    /**
     * @brief Collect values on end nodes into a dense vector.
     * @param[out] data Column vector sized to the total number of end-node values.
     */
    void getEndValues(Eigen::Matrix<T, Eigen::Dynamic, 1> &data);

    /**
     * @brief Set end-node values from a dense vector.
     * @param[in] data Column vector holding values; its size must match.
     */
    void setEndValues(Eigen::Matrix<T, Eigen::Dynamic, 1> &data);

    /**
     * @brief Persist the tree (binary).
     * @param file Output filename.
     */
    void saveTree(const std::string &file);

    /**
     * @brief Persist the tree (text).
     * @param file Output filename.
     */
    void saveTreeTXT(const std::string &file);

    /**
     * @brief Load the tree (binary).
     * @param file Input filename.
     */
    void loadTree(const std::string &file);

    /**
     * @brief Load the tree (text).
     * @param file Input filename.
     */
    void loadTreeTXT(const std::string &file);

    /** @name In-place algebra on the represented function */
    ///@{

    /// Square the function pointwise.
    void square();

    /// Raise the function to power p pointwise.
    void power(double p);

    /// Multiply the function by a scalar c.
    void rescale(T c);

    /// Normalize the function to unit norm (when meaningful).
    void normalize();

    /// Compute this := this + c * inp (alloc/refine as needed).
    void add(T c, FunctionTree<D, T> &inp);

    /// In-place addition on the existing grid only.
    void add_inplace(T c, FunctionTree<D, T> &inp);

    /// Compute this := this + c * |inp| (absolute values).
    void absadd(T c, FunctionTree<D, T> &inp);

    /// Compute this := this * (c * inp) pointwise.
    void multiply(T c, FunctionTree<D, T> &inp);

    /// Apply a scalar-to-scalar map pointwise.
    void map(FMap<T, T> fmap);

    ///@}

    /**
     * @brief Number of memory chunks reserved for nodes.
     * @return Total chunk count.
     */
    int getNChunks() { return this->getNodeAllocator().getNChunks(); }

    /**
     * @brief Number of memory chunks currently in use.
     * @return Used chunk count.
     */
    int getNChunksUsed() { return this->getNodeAllocator().getNChunksUsed(); }

    /**
     * @brief Prune small contributions and optionally refine slightly.
     * @param prec Threshold used for pruning.
     * @param splitFac Optional split factor for balancing.
     * @param absPrec If true, use absolute thresholding.
     * @return Number of nodes removed or affected.
     */
    int crop(double prec, double splitFac = 1.0, bool absPrec = true);

    /** @name Typed access to nodes */
    ///@{

    /// Get i-th end node cast to FunctionNode (non-const).
    FunctionNode<D, T> &getEndFuncNode(int i) { return static_cast<FunctionNode<D, T> &>(this->getEndMWNode(i)); }

    /// Get i-th root node cast to FunctionNode (non-const).
    FunctionNode<D, T> &getRootFuncNode(int i) { return static_cast<FunctionNode<D, T> &>(this->rootBox.getNode(i)); }

    /// Allocator for generated nodes (non-const).
    NodeAllocator<D, T> &getGenNodeAllocator() { return *this->genNodeAllocator_p; }

    /// Allocator for generated nodes (const).
    const NodeAllocator<D, T> &getGenNodeAllocator() const { return *this->genNodeAllocator_p; }

    /// Get i-th end node cast to FunctionNode (const).
    const FunctionNode<D, T> &getEndFuncNode(int i) const { return static_cast<const FunctionNode<D, T> &>(this->getEndMWNode(i)); }

    /// Get i-th root node cast to FunctionNode (const).
    const FunctionNode<D, T> &getRootFuncNode(int i) const { return static_cast<const FunctionNode<D, T> &>(this->rootBox.getNode(i)); }

    ///@}

    /**
     * @brief Delete nodes that were generated during the last build/refine step.
     * @details Restores the tree to the pre-generation state without touching
     * persisted nodes and data.
     */
    void deleteGenerated();

    /**
     * @brief Delete generated nodes and their generated parents if they became empty.
     */
    void deleteGeneratedParents();

    /**
     * @brief Build a flat view of the coefficient storage.
     *
     * @param[out] coefs      Pointers to coefficient blocks per node.
     * @param[out] indices    Node indices mapped to a compact integer id.
     * @param[out] parent_indices Parent ids matching indices.
     * @param[out] scalefac   Per-node scale factors (e.g. for normalization).
     * @param[out] max_index  Maximum assigned compact id.
     * @param[in]  refTree    Reference tree defining traversal order.
     * @param[in]  refNodes   Optional explicit node list to follow.
     *
     * @details Intended for exporting the tree into custom linear algebra
     * back-ends or checkpoint formats.
     */
    void makeCoeffVector(std::vector<T *> &coefs,
                         std::vector<int> &indices,
                         std::vector<int> &parent_indices,
                         std::vector<double> &scalefac,
                         int &max_index,
                         MWTree<D, double> &refTree,
                         std::vector<MWNode<D, double> *> *refNodes = nullptr);

    /**
     * @brief Reconstruct a tree topology from a coefficient vector.
     * @param refTree  Reference topology to follow.
     * @param coefpVec Pointers to coefficient blocks.
     * @param ix2coef  Mapping from node compact id to coefpVec index.
     * @param absPrec  Threshold for adaptive creation.
     * @param mode     Creation mode: "adaptive" or fixed variants.
     */
    void makeTreefromCoeff(MWTree<D, double> &refTree,
                           std::vector<T *> coefpVec,
                           std::map<int, int> &ix2coef,
                           double absPrec,
                           const std::string &mode = "adaptive");

    /**
     * @brief Append topology from another tree (no coefficients copied).
     * @param inTree Input tree.
     */
    void appendTreeNoCoeff(MWTree<D, double> &inTree);

    /// @overload
    void appendTreeNoCoeff(MWTree<D, ComplexDouble> &inTree);

    /**
     * @brief Copy topology and coefficients from a real-valued tree.
     * @param inTree Source tree.
     */
    void CopyTree(FunctionTree<D, double> &inTree);

    /**
     * @brief Move all node coefficients to a bank and remove them from nodes.
     * @return Number of nodes affected.
     */
    int saveNodesAndRmCoeff();

    /**
     * @brief Deep-copy entire tree into out (topology and data).
     * @param out Destination tree pointer (must be non-null and compatible).
     */
    void deep_copy(FunctionTree<D, T> *out);

    /**
     * @brief Extract real part into a newly allocated real tree.
     * @return Pointer to a new FunctionTree of type double.
     */
    FunctionTree<D, double> *Real();

    /**
     * @brief Extract imaginary part into a newly allocated real tree.
     * @return Pointer to a new FunctionTree of type double.
     */
    FunctionTree<D, double> *Imag();

    /** @name Real/complex conversion helpers */
    ///@{
    void CopyTreeToComplex(FunctionTree<3, ComplexDouble> *&out);
    void CopyTreeToComplex(FunctionTree<2, ComplexDouble> *&out);
    void CopyTreeToComplex(FunctionTree<1, ComplexDouble> *&out);
    void CopyTreeToReal(FunctionTree<3, double> *&out); // for testing
    ///@}

protected:
    /// Allocator for generated nodes.
    std::unique_ptr<NodeAllocator<D, T>> genNodeAllocator_p{nullptr};

    /// Print a short, human-readable description of the tree.
    std::ostream &print(std::ostream &o) const override;

    /// Allocate and initialize root nodes according to the MRA.
    void allocRootNodes();
};

} // namespace mrcpp