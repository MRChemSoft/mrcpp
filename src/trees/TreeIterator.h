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

/**
 * @file TreeIterator.h
 * @brief Iteration helpers for traversing multiwavelet trees.
 *
 * @details
 * This header provides a depth-aware iterator over the nodes of a @ref MWTree.
 * It supports different **traversal directions** and **node-ordering schemes**,
 * selected via constants defined in @c MRCPP/constants.h:
 * - Traversal mode: @c TopDown or @c BottomUp
 *   In the @c TopDown mode, one iterates from the first root node and recursively 
 *   over the children
 *   In the @c BottomUp mode, one first traverses the tree all the way down to the 
 *   leaves and then starts iteratig from there
 * - Iterator type:  @c Lebesgue (Z-order) or @c Hilbert
 *
 * The iterator yields @ref MWNode instances in the requested sequence determined by 
 * the parameters above
 *
 * The file contains two classes: @ref TreeIterator and @ref IteratorNode.
 * The @ref TreeIterator is the main interface for users, while the @ref IteratorNode 
 * is mainly a placeholder for a few node-specific flags.
 *
 * @par Example
 * @code{.cpp}
 * using namespace mrcpp;
 * TreeIterator<3,double> it(tree, TopDown, Lebesgue);
 * it.setReturnGenNodes(true);   // include generated nodes
 * it.setMaxDepth(5);            // restrict to depth <= 5
 *
 * while (it.next()) {
 *     MWNode<3,double> &nd = it.getNode();
 *     // ... inspect nd, read coefficients/norms, etc.
 * }
 * @endcode
 */

#pragma once

#include "MRCPP/constants.h"
#include "MRCPP/mrcpp_declarations.h"

namespace mrcpp {

/**
 * @class TreeIterator
 *
 * @brief Iterator for traversing an @ref MWTree.
 *
 * @tparam D Spatial dimension (1, 2,  or 3)
 * @tparam T Coefficient type (e.g. double,  ComplexDouble) 
 *
 * @details
 * The iterator traverses the tree starting the root node(s), producing nodes
 * according to:
 * - a **traversal direction** ( @c TopDown or @c BottomUp), and
 * - an **ordering scheme** within siblings ( @c Lebesgue or @c Hilbert).
 *
 * The behavior can be refined with:
 * - @ref setReturnGenNodes() to toggle inclusion of generated (non-leaf) nodes,
 * - @ref setMaxDepth() to limit the traversal depth,
 * - @ref setTraverse() / @ref setIterator() to change policies at runtime.
 *
 * The iteration state is represented by an internal linked stack of
 * @ref IteratorNode instances.
 */
template <int D, typename T> class TreeIterator {
public:
    /**
     * @brief Construct a detached iterator (no tree bound yet).
     *
     * @param traverse Traversal mode (e.g., @c TopDown or @c BottomUp).
     * @param iterator Node-ordering mode (e.g., @c Lebesgue or @c Hilbert).
     *
     * @note Call @ref init() before the first @ref next() if you use this ctor.
     */
    TreeIterator(int traverse = TopDown, int iterator = Lebesgue);
    /**
     * @brief Construct an iterator bound to a tree.
     *
     * @param tree      Tree to traverse.
     * @param traverse  Traversal mode (e.g., @c TopDown or @c BottomUp).
     * @param iterator  Node-ordering mode (e.g., @c Lebesgue or @c Hilbert).
     */
    TreeIterator(MWTree<D, T> &tree, int traverse = TopDown, int iterator = Lebesgue);
    /// @brief Destructor (releases internal traversal state).
    virtual ~TreeIterator();
    void setReturnGenNodes(bool i = true) { this->returnGenNodes = i; } ///< @param i If true, generated nodes are included in the sequence.
    void setMaxDepth(int depth) { this->maxDepth = depth; } ///< @param depth Non-negative maximum depth; if negative, no limit is applied.
    void setTraverse(int traverse);///< @param traverse set Traversal mode (@c TopDown or @c BottomUp).
    void setIterator(int iterator);///< @param iterator set Iterator type (@c Lebesgue or @c Hilbert).
    MWNode<D, T> &getNode() { return *this->state->node; } ///< @return Reference to the node yielded by the last successful @ref next() / @ref nextParent().
    /**
     * @brief Bind the iterator to a tree and reset traversal state.
     *
     * @param tree Tree to traverse.
     */
    void init(MWTree<D, T> &tree);
    /**
     * @brief Advance to the next node according to the current policy.
     *
     * @return @c true if a node is available (use @ref getNode()), @c false when finished.
     *
     * @details
     * if the current @ref IteratorNode is null, return false.
     * In @c TopDown mode, try to return the current node first.
     * If successful, return true.
     * If not, check if the current node has children, and try to return
     * the next child node according to the ordering scheme.
     * If successful, return true.
     * If not, try to move to the next root node, and return its first node
     * according to the ordering scheme.
     * If successful, return true.
     * If not, in @c BottomUp mode, try to return the current node.
     * If successful, return true.
     * If not, remove the current state and recur invoking a new @ref next().
     */
    bool next();
    /**
     * @brief Advance to the next parent node according to the current policy.
     *
     * @return @c true if the parent node is available, @c false when finished.
     *
     * @details
     * Returns the current node or the parent of the current node. The logic makes sure the correct
     * parent is returned according to the traversal mode and ordering scheme. In case of PBC calculations,
     * the parent may be above the root nodes defining the unit cell.
     */
    bool nextParent();

    friend class IteratorNode<D, T>;

protected:
    int root{0};                ///< Index of the current root box.
    int nRoots{0};              ///< Number of root boxes in the tree.
    int mode{TopDown};          ///< Traversal mode (@c TopDown or @c BottomUp).
    int type{Lebesgue};         ///< Iterator type (@c Lebesgue or @c Hilbert).
    int maxDepth{-1};           ///< Max depth limit; negative means unlimited.
    bool returnGenNodes{true};  ///< If @c true, also return generated (non-leaf) nodes.
    IteratorNode<D, T> *state{nullptr};         ///< Current traversal frame.
    IteratorNode<D, T> *initialState{nullptr};  ///< Initial frame for the current root.

    int getChildIndex(int i) const; ///< @brief Map logical child order [0..2^D) to actual child index based on @ref type.
/**
 * @name try... methods
 * @brief The following methods test if the node of a given type should be returned.
 * @details In addition to returning @c true or @c false, these methods also update the internal
 * traversal state accordingly.
 * @{
 */
    bool tryParent(); ///< @return @c true if the parent node should be returned.
    bool tryChild(int i);///< @return @c true if the child at index @p i should be returned.
    bool tryNode(); ///< @return @c true if the current node shuld be returned.
    bool tryNextRoot(); ///< @return @c true if the next root node should be returned.
    bool tryNextRootParent(); ///< @return @c true if the parent of the next root node is available and should be returned.
/** @} */
    void removeState(); ///< @brief Remove the current traversal frame from the stack.
    bool checkDepth(const MWNode<D, T> &node) const; ///< @return @c true if the node is within the max depth limit.
    bool checkGenerated(const MWNode<D, T> &node) const; ///< @return @c true if the generated nodes should be included.
};

/**
 * @class IteratorNode
 * @brief Iterator representing a node in the traversal stack.
 *
 * @tparam D Spatial dimension (1, 2,  or 3)
 * @tparam T Coefficient type (e.g. double,  ComplexDouble) 
 *
 * @details
 * This is an internal placeholder which contains both the pointer to the actual node to return and 
 * flags to determine if itself, its parent and its children have been already returned. 
  * It contains:
 * - a pointer to the node,
 * - a link to the next node in the stack
 * - completion flags for the current node, its parent, and its children.
 */
template <int D, typename T> class IteratorNode final {
public:
    MWNode<D, T> *node;              ///< Current node.
    IteratorNode<D, T> *next;        ///< Next node in the stack.
    bool doneNode;                   ///< Whether the node itself has been used.
    bool doneParent;                 ///< Whether the parent node has been used.
    bool doneChild[1 << D];          ///< Whether each child has been used.

    /**
     * @brief Construct a new iterator
     *
     * @param nd Pointer to the MW node represented by this frame.
     * @param nx Link to the next iterator (can be @c nullptr).
     */
    IteratorNode(MWNode<D, T> *nd, IteratorNode<D, T> *nx = nullptr);

    /// @brief Recursively delete the linked iterators that follow this one.
    ~IteratorNode() { delete this->next; }
};

} // namespace mrcpp