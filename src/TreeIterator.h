#pragma once

#include "constants.h"
#include "mrcpp_declarations.h"

template<int D> class IteratorNode;

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
    bool returnGenNodes;
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
class IteratorNode {
public:
    MWNode<D> *node;
    IteratorNode<D> *next;
    bool doneNode;
    bool doneChild[1 << D];

    IteratorNode(MWNode<D> *nd, IteratorNode<D> *nx = 0);
    ~IteratorNode() { delete this->next; }
};

