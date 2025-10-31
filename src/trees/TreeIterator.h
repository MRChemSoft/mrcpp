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
 * This header provides a generic depth-aware iterator over @ref MWTree nodes.
 * It supports different **traversal directions** and **node-ordering schemes**,
 * selected via constants defined in @c MRCPP/constants.h:
 * - Traversal mode: @c TopDown or @c BottomUp
 * - Iterator type:  @c Lebesgue (Z-order) or @c Hilbert (space-filling)
 *
 * The iterator yields @ref MWNode instances from one or more root nodes,
 * honoring a user-provided maximum depth and whether *generated* (non-end)
 * nodes should be returned.
 *
 * @par Example
 * @code{.cpp}
 * using namespace mrcpp;
 * TreeIterator<3,double> it(tree, TopDown, Lebesgue);
 * it.setReturnGenNodes(true);   // include generated/branch nodes
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
 * @brief Stateful iterator for traversing an @ref MWTree.
 *
 * @tparam D Spatial dimensionality (1, 2, or 3).
 * @tparam T Coefficient type (e.g., @c double or @c ComplexDouble).
 *
 * @details
 * The iterator walks the tree starting from each root node, producing nodes
 * according to:
 * - a **traversal direction** (@c TopDown or @c BottomUp), and
 * - an **ordering scheme** within siblings (@c Lebesgue or @c Hilbert).
 *
 * The behavior can be refined with:
 * - @ref setReturnGenNodes() to toggle inclusion of generated (non-leaf) nodes,
 * - @ref setMaxDepth() to limit the traversal depth,
 * - @ref setTraverse() / @ref setIterator() to change policies at runtime.
 *
 * The iteration state is represented by a small internal linked stack of
 * @ref IteratorNode frames.
 */
template <int D, typename T> class TreeIterator {
public:
    /**
     * @brief Construct a detached iterator (no tree bound yet).
     * @param traverse Traversal mode (e.g., @c TopDown or @c BottomUp).
     * @param iterator Node-ordering mode (e.g., @c Lebesgue or @c Hilbert).
     *
     * @note Call @ref init() before the first @ref next() if you use this ctor.
     */
    TreeIterator(int traverse = TopDown, int iterator = Lebesgue);

    /**
     * @brief Construct an iterator bound to a tree.
     * @param tree      Tree to traverse.
     * @param traverse  Traversal mode (e.g., @c TopDown or @c BottomUp).
     * @param iterator  Node-ordering mode (e.g., @c Lebesgue or @c Hilbert).
     */
    TreeIterator(MWTree<D, T> &tree, int traverse = TopDown, int iterator = Lebesgue);

    /// @brief Destructor (releases internal traversal state).
    virtual ~TreeIterator();

    /**
     * @brief Include/exclude generated (non-end) nodes in the iteration stream.
     * @param i If @c true, generated nodes are returned by @ref next().
     *          If @c false, only end (leaf) nodes are produced.
     */
    void setReturnGenNodes(bool i = true) { this->returnGenNodes = i; }

    /**
     * @brief Set maximum depth measured from the root scale.
     * @param depth Non-negative maximum depth; if negative, no limit is applied.
     */
    void setMaxDepth(int depth) { this->maxDepth = depth; }

    /**
     * @brief Change traversal mode at runtime.
     * @param traverse @c TopDown or @c BottomUp (see @c MRCPP/constants.h).
     * @warning Changing mode invalidates in-flight assumptions; call before @ref init().
     */
    void setTraverse(int traverse);

    /**
     * @brief Change sibling-ordering policy at runtime.
     * @param iterator @c Lebesgue or @c Hilbert (see @c MRCPP/constants.h).
     * @warning Changing mode invalidates in-flight assumptions; call before @ref init().
     */
    void setIterator(int iterator);

    /**
     * @brief Bind the iterator to a tree and reset traversal state.
     * @param tree Tree to traverse.
     */
    void init(MWTree<D, T> &tree);

    /**
     * @brief Advance to the next node according to the current policy.
     * @return @c true if a node is available (use @ref getNode()), @c false when finished.
     */
    bool next();

    /**
     * @brief Move the cursor to the parent of the current node (if any).
     * @return @c true if the parent exists and becomes current, otherwise @c false.
     */
    bool nextParent();

    /**
     * @brief Access the current node.
     * @return Reference to the node yielded by the last successful @ref next() / @ref nextParent().
     */
    MWNode<D, T> &getNode() { return *this->state->node; }

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

    /// @brief Map logical child order [0..2^D) to physical child index based on @ref type.
    int getChildIndex(int i) const;

    /// @name Traversal helpers
    ///@{
    bool tryParent();
    bool tryChild(int i);
    bool tryNode();
    bool tryNextRoot();
    bool tryNextRootParent();
    void removeState();
    bool checkDepth(const MWNode<D, T> &node) const;
    bool checkGenerated(const MWNode<D, T> &node) const;
    ///@}
};

/**
 * @class IteratorNode
 * @brief Lightweight frame holding traversal state for one MW node.
 *
 * @tparam D Spatial dimensionality (1, 2, or 3).
 * @tparam T Coefficient type (e.g., @c double or @c ComplexDouble).
 *
 * @details
 * The iterator maintains a small linked list (stack) of these frames while
 * walking the tree. Each frame keeps:
 * - a pointer to the node,
 * - a link to the previous frame,
 * - completion flags for the current node, its parent, and its children.
 */
template <int D, typename T> class IteratorNode final {
public:
    MWNode<D, T> *node;              ///< Current node.
    IteratorNode<D, T> *next;        ///< Previous frame in the stack.
    bool doneNode;                   ///< Whether the node itself has been yielded.
    bool doneParent;                 ///< Whether the parent transition has been attempted.
    bool doneChild[1 << D];          ///< Whether each child has been attempted.

    /**
     * @brief Construct a traversal frame.
     * @param nd Pointer to the MW node represented by this frame.
     * @param nx Link to the previous frame (can be @c nullptr).
     */
    IteratorNode(MWNode<D, T> *nd, IteratorNode<D, T> *nx = nullptr);

    /// @brief Recursively delete the linked frames that follow this one.
    ~IteratorNode() { delete this->next; }
};

} // namespace mrcpp