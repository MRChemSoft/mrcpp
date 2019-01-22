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

#include "constants.h"
#include "mrcpp_declarations.h"

namespace mrcpp {

template<int D>
class TreeIterator {
public:
    TreeIterator(int dir = TopDown);
    virtual ~TreeIterator();

    void setReturnGenNodes(bool i = true) { this->returnGenNodes = i; }
    void setMaxDepth(int depth) { this->maxDepth = depth; }

    bool next();
    MWNode<D> &getNode() { return *this->state->node; }

    friend class IteratorNode<D>;

protected:
    int root;
    int nRoots;
    int mode;
    int maxDepth;
    bool returnGenNodes{true};
    IteratorNode<D> *state;
    IteratorNode<D> *initialState;

    virtual int getChildIndex(int i) const = 0;

    void init(MWTree<D> *tree);
    bool tryNode();
    bool tryChild(int i);
    bool tryNextRoot();
    void removeState();
    void setDirection(int dir);
    bool checkDepth(const MWNode<D> &node) const;
    bool checkGenerated(const MWNode<D> &node) const;
};

template<int D>
class IteratorNode final {
public:
    MWNode<D> *node;
    IteratorNode<D> *next;
    bool doneNode;
    bool doneChild[1 << D];

    IteratorNode(MWNode<D> *nd, IteratorNode<D> *nx = 0);
    ~IteratorNode() { delete this->next; }
};

}
