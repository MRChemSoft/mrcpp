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

template <int D> class TreeIterator {
public:
    TreeIterator(int traverse = TopDown, int iterator = Lebesgue);
    TreeIterator(MWTree<D> &tree, int traverse = TopDown, int iterator = Lebesgue);
    virtual ~TreeIterator();

    void setReturnGenNodes(bool i = true) { this->returnGenNodes = i; }
    void setMaxDepth(int depth) { this->maxDepth = depth; }
    void setTraverse(int traverse);
    void setIterator(int iterator);

    void init(MWTree<D> &tree);
    bool next();
    bool nextParent();
    MWNode<D> &getNode() { return *this->state->node; }

    friend class IteratorNode<D>;

protected:
    int root;
    int nRoots;
    int mode;
    int type;
    int maxDepth;
    bool returnGenNodes{true};
    IteratorNode<D> *state;
    IteratorNode<D> *initialState;

    int getChildIndex(int i) const;

    bool tryParent();
    bool tryChild(int i);
    bool tryNode();
    bool tryNextRoot();
    bool tryNextRootParent();
    void removeState();
    bool checkDepth(const MWNode<D> &node) const;
    bool checkGenerated(const MWNode<D> &node) const;
};

template <int D> class IteratorNode final {
public:
    MWNode<D> *node;
    IteratorNode<D> *next;
    bool doneNode;
    bool doneParent;
    bool doneChild[1 << D];

    IteratorNode(MWNode<D> *nd, IteratorNode<D> *nx = nullptr);
    ~IteratorNode() { delete this->next; }
};

} // namespace mrcpp
