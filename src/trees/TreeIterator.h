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

#include "MRCPP/constants.h"
#include "MRCPP/mrcpp_declarations.h"

namespace mrcpp {

template <int D, typename T> class TreeIterator {
public:
    TreeIterator(int traverse = TopDown, int iterator = Lebesgue);
    TreeIterator(MWTree<D, T> &tree, int traverse = TopDown, int iterator = Lebesgue);
    virtual ~TreeIterator();

    void setReturnGenNodes(bool i = true) { this->returnGenNodes = i; }
    void setMaxDepth(int depth) { this->maxDepth = depth; }
    void setTraverse(int traverse);
    void setIterator(int iterator);

    void init(MWTree<D, T> &tree);
    bool next();
    bool nextParent();
    MWNode<D, T> &getNode() { return *this->state->node; }

  friend class IteratorNode<D, T>;

protected:
    int root;
    int nRoots;
    int mode;
    int type;
    int maxDepth;
    bool returnGenNodes{true};
    IteratorNode<D, T> *state;
    IteratorNode<D, T> *initialState;

    int getChildIndex(int i) const;

    bool tryParent();
    bool tryChild(int i);
    bool tryNode();
    bool tryNextRoot();
    bool tryNextRootParent();
    void removeState();
    bool checkDepth(const MWNode<D, T> &node) const;
    bool checkGenerated(const MWNode<D, T> &node) const;
};

template <int D, typename T> class IteratorNode final {
public:
    MWNode<D, T> *node;
    IteratorNode<D, T> *next;
    bool doneNode;
    bool doneParent;
    bool doneChild[1 << D];

    IteratorNode(MWNode<D, T> *nd, IteratorNode<D, T> *nx = nullptr);
    ~IteratorNode() { delete this->next; }
};

} // namespace mrcpp
