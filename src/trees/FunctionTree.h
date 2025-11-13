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
 * @class FunctionTree
 * @tparam D Spatial dimension (1, 2, or 3).
 * @tparam T Coefficient type (e.g. double, ComplexDouble).
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
     * @param mra Multi-resolution analysis (domain and basis function used).
     * @param name Optional textual name of the function.
     *
     * @note This will only create the object. To compute coefficients,
     *       use a tree builder or calculator afterwards (e.g., projectors).
     */
    FunctionTree(const MultiResolutionAnalysis<D> &mra, const std::string &name)
            : FunctionTree(mra, nullptr, name) {}

    /**
     * @brief Construct a tree bound to an MRA with optional shared memory and name.
     * @param[in] mra: Which MRA the function is defined
     * @param[in] sh_mem: Pointer to MPI shared memory block
     * @param name Optional textual name of the function
     * 
     * @returns New FunctionTree object
     *
     * @details Constructs an uninitialized tree, containing only empty root nodes. 
     * If a shared memory pointer is provided the tree will be allocated in this
     * shared memory window, otherwise it will be local to each MPI process.
     */
    FunctionTree(const MultiResolutionAnalysis<D> &mra,
                 SharedMemory<T> *sh_mem = nullptr,
                 const std::string &name = "nn");

    
    FunctionTree(const FunctionTree<D, T> &tree) = delete;
    FunctionTree<D, T> &operator=(const FunctionTree<D, T> &tree) = delete;

    /// FunctionTree destructor
    ~FunctionTree() override;

    /**
     * @brief Integrate the represented function over the MRA box
     * @returns Integral of the function over the entire computational domain 
     */
    T integrate() const;

    /**
     * @brief Integrate a representable function using the this tree's grid
     * @param[in] f RepresentableFunction used as integrand partner
     * @returns Integral of the representable function
     *
     * @details You can evaluate the integral of any representable function
     * over the most refined scale of the 'this' FunctionTree's grid.
     */
    double integrateEndNodes(RepresentableFunction_M &f);

    /**
     * @brief Evaluate with high accuracy at a given coordinate.
     * @param r Physical coordinate.
     * @return Function value.
     *
     * @details May be more expensive than evalf due to stricter handling.
     */
    

    /**
     * @brief Fast but approximate evaluation of this FunctionTree at a given coordinate.
     * @param[in] r: Cartesian coordinate to be evaluated
     * @return Approximate Function value
     * 
     * @note This will only evaluate the _scaling_ part of the
     *       leaf nodes in the tree, which means that the function
     *       values will not be fully accurate.
     *       This is done to allow a fast and const function evaluation
     *       that can be done in OMP parallel. If you want to include
     *       also the _final_ wavelet part you can call the corresponding
     *       evalf_precise function, _or_ you can manually extend
     *       the MW grid by one level before evaluating, using
     *       `mrcpp::refine_grid(tree, 1)`
     */
    T evalf(const Coord<D> &r) const override;

    /** 
     * @brief Slow but high-accuracy , evaluation of the function at a given coordinate
     * @param[in] r: Cartesian coordinate to be evaluated
     * @returns Exact value of this FunctionTree in the point r
     * @note This will evaluate the _true_ value (scaling + wavelet) of the
     *       leaf nodes in the tree. This requires an on-the-fly MW transform
     *       on the node which makes this function slow and non-const. If you
     *       need fast evaluation, use refine_grid(tree, 1) first, and then
     *       evalf.     
     */
    T evalf_precise(const Coord<D> &r);

    /**
     * @brief Number of generated (non-root) nodes currently alive
     * @return Count of nodes (managed by the generated-node allocator)
     */
    int getNGenNodes() const { return getGenNodeAllocator().getNNodes(); }

    /**
     * @brief Collect values on end nodes into a dense Eigen type vector
     * @param[out] data Column vector sized to the total number of end-node values
     */
    void getEndValues(Eigen::Matrix<T, Eigen::Dynamic, 1> &data);

    /**
     * @brief Set end-node values as the components of the dense Eigen type vector
     * @param[in] data Column vector holding values; its size must match
     */
    void setEndValues(Eigen::Matrix<T, Eigen::Dynamic, 1> &data);

    /**
     * @brief Write the tree to disk in text/ASCII format in a representation
     *   using MADNESS conventions for n, l and index order.
     * @param file Output filename.
     */
    void saveTreeTXT(const std::string &file);

    /**
     * @brief Write the tree structure to disk, for later use
     * @param[in] file: File name, will get ".tree" extension
     */
    void saveTree(const std::string &file);

    /**
     * @brief Read a previously stored tree structure from disk
     * @param[in] file File name, will get ".tree" extension
     * @note This tree must have the exact same MRA the one that was saved
     */
    void loadTree(const std::string &file);

    /**
     * @brief Read a previously stored tree assuming text/ASCII format, using MADNESS conventions (n, l and index order)
     * @param[in] file Input filename
     * @note Make sure that the MRA of this tree matches the one used to create the file
     */
    void loadTreeTXT(const std::string &file);


    /** @brief In-place square of MW function representations, fixed grid
     *
     * @details The leaf node point values of the function will be in-place
     * squared, no grid refinement.
     *
     */
    void square(); 
    /// Raise the function to power p pointwise.
    /**
     * @brief In-place power of MW function representations, fixed grid
     * 
     * @param p Exponent
     * 
     * @details The leaf node point values of the function will be in-place raised
     * to the given power, no grid refinement.
     */
    void power(double p);
    /**
     * @brief In-place multiplication by a scalar, fixed grid
     * 
     * @param c Scaling factor (with the same data type as the coefficients)
     * 
     * @details The leaf node point values of the function will be
     * in-place multiplied by the given coefficient, no grid refinement.
     */
    void rescale(T c);  
    void normalize();  ///< In-place rescaling by a function norm \f$ ||f||^{-1} \f$, fixed grid
    
    /**
     * @brief this + inp (fixed grid)
     * 
     * @param c: Numerical coefficient of input function
     * @param[in] inp: Input function to be added on this FunctionTree
     * 
     * @details The input function will be added in-place on the current grid of
     * the function, i.e. no further grid refinement. Addition done within the MW representations.
     */
    void add(T c, FunctionTree<D, T> &inp);

    /**
     * @brief this + inp (uniting the two grids)
     
     * @param c: Numerical coefficient of input function
     * @param[in] inp: Input function to be added on this FunctionTree
     * 
     * @details The input function will be added to the union of the current grid of
     * and input the function grid. Addition done within the MW representations.
     */
    void add_inplace(T c, FunctionTree<D, T> &inp);

    /**
     * @brief this + abs(inp) (fixed grid)
     
     * @param c: Numerical coefficient of input function
     * @param[in] inp: Input function to be added on this FunctionTree
     * 
     * @details The absolute value of input function will be added in-place on the current grid of the output
     * function, i.e. no further grid refinement. Addition done within the MW representations.
     */
    void absadd(T c, FunctionTree<D, T> &inp);

    /**
     * @brief this * (c * inp), fixed grid
     * @param c: Numerical coefficient of input function
     * @param[in] inp: Input function to be multiplied with this FunctionTree
     * 
     * @details The input function will be multiplied in-place on the current grid
     * of the function, i.e. no further grid refinement.
     */
    void multiply(T c, FunctionTree<D, T> &inp);

    /** 
     * @brief In-place mapping with a predefined function f(x), fixed grid
     * 
     * @param[in] fmap: mapping function
     *
     * @details The input function will be mapped in-place on the current grid
     * of the function, i.e. no further grid refinement.
     */
    void map(FMap<T, T> fmap);




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

    /** @brief Reduce the precision of the tree by deleting nodes
     *
     * @param prec: New precision criterion
     * @param splitFac: Splitting factor: 1, 2 or 3
     * @param absPrec: Use absolute precision
     *
     * @details This will run the tree building algorithm in "reverse", starting
     * from the leaf nodes, and perform split checks on each node based on the given
     * precision and the local wavelet norm.
     *
     * @note The splitting factor appears in the threshold for the wavelet norm as
     * \f$ ||w|| < 2^{-sn/2} ||f|| \epsilon \f$. In principal, `s` should be equal
     * to the dimension; in practice, it is set to `s=1`.
     */
    int crop(double prec, double splitFac = 1.0, bool absPrec = true);

    /** @name Typed access to nodes */
    ///@{

    /// @return i-th end node cast to FunctionNode 
    FunctionNode<D, T> &getEndFuncNode(int i) { return static_cast<FunctionNode<D, T> &>(this->getEndMWNode(i)); }

    /// @return i-th root node cast to FunctionNode
    FunctionNode<D, T> &getRootFuncNode(int i) { return static_cast<FunctionNode<D, T> &>(this->rootBox.getNode(i)); }

    /// @return Allocator for generated nodes 
    NodeAllocator<D, T> &getGenNodeAllocator() { return *this->genNodeAllocator_p; }

    /// @return Allocator for generated nodes ìì
    const NodeAllocator<D, T> &getGenNodeAllocator() const { return *this->genNodeAllocator_p; }

    /// @return i-th end node cast to FunctionNode 
    const FunctionNode<D, T> &getEndFuncNode(int i) const { return static_cast<const FunctionNode<D, T> &>(this->getEndMWNode(i)); }

    /// @return i-th root node cast to FunctionNode
    const FunctionNode<D, T> &getRootFuncNode(int i) const { return static_cast<const FunctionNode<D, T> &>(this->rootBox.getNode(i)); }

    ///@}

    /**
     * @brief Delete nodes that were generated during the last build/refine step.
     * @details Restores the tree to the pre-generation state without touching
     * persisted nodes and data, lowering the resolution to the previous scale
     */
    void deleteGenerated();

    /**
     * @brief Delete generated nodes and their generated parents if they became empty.
     */
    void deleteGeneratedParents();

    /**
     * @brief Will fill the first 4 vectors with the coefficient pointers, indices, parent indices and scale factors
     *
     * @param[out] coefs      Pointers to coefficient blocks per node.
     * @param[out] indices    Node indices mapped to a compact integer id.
     * @param[out] parent_indices Parent ids matching indices.
     * @param[out] scalefac   Per-node scale factors (e.g. for normalization).
     * @param[out] max_index  Maximum assigned compact id.
     * @param[in]  refTree    Reference tree defining traversal order.
     * @param[in]  refNodes   Optional explicit node list to follow.
     *
     * @details Traverse tree using BFS and returns an array with the address of the coefs.
     * Also returns an array with the corresponding indices defined as the
     * values of serialIx in refTree, and an array with the indices of the parent.
     * Set index -1 for nodes that are not present in refTree
     * Set parent_indices as -2 for nodes that are not present in refTree
     * Set scalefac as 1.0 for nodes that are not present in refTree
     * Intended for exporting the tree into custom linear algebra
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
     * @param[out] refTree  Reference topology to follow.
     * @param coefpVec Pointers to coefficient blocks.
     * @param ix2coef  Mapping from node compact id to coefpVec index.
     * @param absPrec  Threshold for adaptive creation.
     * @param mode     Creation mode: "adaptive" or fixed variants.
     * 
     * @details Traverse tree using DFS and reconstruct it using node info from the
     * reference tree and a list of coefficients.
     * It is the reference tree (refTree) which is traversed, but one does not descend
     * into children if the norm of the tree is smaller than absPrec.
     */
    void makeTreefromCoeff(MWTree<D, double> &refTree,
                           std::vector<T *> coefpVec,
                           std::map<int, int> &ix2coef,
                           double absPrec,
                           const std::string &mode = "adaptive");


    /**
     * @brief Append topology from another tree with real-type coefficients (no coefficients copied)
     * @param[in] inTree Input tree.
     * 
     * @note It will append only the nodes structure, without copying any coefficient data,
     * therefore it won't matter the datatype of the input tree for the result.
     */
    void appendTreeNoCoeff(MWTree<D, double> &inTree);

    /**
     * @brief Append topology from another tree with real-type coefficients (no coefficients copied)
     * @param[in] inTree Input tree.
     * 
     * @note It will append only the nodes structure, without copying any coefficient data,
     * therefore it won't matter the datatype of the input tree for the result.
     */
    void appendTreeNoCoeff(MWTree<D, ComplexDouble> &inTree);

    /**
     * @brief Copy topology AND coefficients from a real-valued tree
     * @param[in] inTree Source tree.
     * 
     * @note The copy process is a shallow copy for the coefficients, i.e.
     * the new tree nodes will point to the same coefficient blocks as the input tree.
     * Therefore, modifying the coefficients in one tree will affect the other.
     */
    void CopyTree(FunctionTree<D, double> &inTree);

    /**
     * @brief Move all node coefficients to a bank and remove them from nodes
     * @return Number of nodes affected.
     */
    int saveNodesAndRmCoeff();

    /**
     * @brief Deep-copy entire tree into out (topology and data)
     * @param[out] out Destination tree pointer (must be non-null and compatible)
     * 
     * @details Exact copy without any binding between old and new tree
     */
    void deep_copy(FunctionTree<D, T> *out);

    /**
     * @brief Extract real part into a newly allocated real tree
     * @return Pointer to a new FunctionTree of type double
     */
    FunctionTree<D, double> *Real();

    /**
     * @brief Extract imaginary part into a newly allocated real tree
     * @return Pointer to a new FunctionTree of type double
     */
    FunctionTree<D, double> *Imag();

    /** @name Real/complex conversion helpers */
    ///@{
    /**
     * @brief Deep-copy this tree into a complex-valued tree.
     * 
     * @param[out] out Destination tree pointer (must be non-null).
     * @details Exact copy into a complex tree, with imaginary parts set to zero
     */
    void CopyTreeToComplex(FunctionTree<3, ComplexDouble> *&out);
    /**
     * @brief Deep-copy this tree into a complex-valued tree.
     * 
     * @param[out] out Destination tree pointer (must be non-null).
     * @details Exact copy into a complex tree, with imaginary parts set to zero
     */
    void CopyTreeToComplex(FunctionTree<2, ComplexDouble> *&out);
    /**
     * @brief Deep-copy this tree into a complex-valued tree.
     * 
     * @param[out] out Destination tree pointer (must be non-null).
     * @details Exact copy into a complex tree, with imaginary parts set to zero
     */
    void CopyTreeToComplex(FunctionTree<1, ComplexDouble> *&out);

    /**
     * @brief Deep-copy this tree into a real-valued tree.
     * 
     * @param[out] out Destination tree pointer (must be non-null).
     * @details Exact copy into a real tree, taking only the real parts
     */
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