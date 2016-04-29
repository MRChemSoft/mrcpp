#ifndef OPERATORNODE_H
#define OPERATORNODE_H

#include "MWNode.h"
#include "OperatorTree.h"

class OperatorNode : public MWNode<2> {
public:
    OperatorTree &getOperTree() { return static_cast<OperatorTree &>(*this->tree); }
    OperatorNode &getOperParent() { return static_cast<OperatorNode &>(*this->parent); }
    OperatorNode &getOperChild(int i) { return static_cast<OperatorNode &>(*this->children[i]); }

    const OperatorTree &getOperTree() const { return static_cast<const OperatorTree &>(*this->tree); }
    const OperatorNode &getOperParent() const { return static_cast<const OperatorNode &>(*this->parent); }
    const OperatorNode &getOperChild(int i) const { return static_cast<const OperatorNode &>(*this->children[i]); }

    friend class MWNode<2>;
    friend class OperatorTree;

protected:
    OperatorNode(OperatorTree &t, const NodeIndex<2> &n);
    OperatorNode(OperatorNode &p, int c);
    OperatorNode(const OperatorNode &n) : MWNode<2>(n) { NOT_IMPLEMENTED_ABORT; }
    OperatorNode& operator=(const OperatorNode &n) { NOT_IMPLEMENTED_ABORT; }
    virtual ~OperatorNode();

    double calcComponentNorm(int i) const;

    void genChildren();
    void createChildren() { MWNode<2>::createChildren(); this->clearIsEndNode(); }
    void deleteChildren() { MWNode<2>::deleteChildren(); this->setIsEndNode(); }

private:
    void createChild(int i);
    void genChild(int i);
};

#endif // OPERATORNODE_H
