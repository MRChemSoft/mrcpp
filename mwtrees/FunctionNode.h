#ifndef FUNCTIONNODE_H
#define FUNCTIONNODE_H

#include "MWNode.h"
#include "FunctionTree.h"

template<int D>
class FunctionNode : public MWNode<D> {
public:
    FunctionTree<D> &getFuncTree() { return static_cast<FunctionTree<D> &>(*this->tree); }
    FunctionNode<D> &getFuncParent() { return static_cast<FunctionNode<D> &>(*this->parent); }
    FunctionNode<D> &getFuncChild(int i) { return static_cast<FunctionNode<D> &>(*this->children[i]); }

    const FunctionTree<D> &getFuncTree() const { return static_cast<const FunctionTree<D> &>(*this->tree); }
    const FunctionNode<D> &getFuncParent() const { return static_cast<const FunctionNode<D> &>(*this->parent); }
    const FunctionNode<D> &getFuncChild(int i) const { return static_cast<const FunctionNode<D> &>(*this->children[i]); }

    void setValues(Eigen::VectorXd &vec);
    void getValues(Eigen::VectorXd &vec);

    friend class FunctionTree<D>;

protected:
    FunctionNode(FunctionTree<D> &t, const NodeIndex<D> &n) : MWNode<D>(t, n) { }
    FunctionNode(FunctionNode<D> &p, int c) : MWNode<D>(p, c) { }
    FunctionNode(const FunctionNode<D> &n) : MWNode<D>(n) { NOT_IMPLEMENTED_ABORT; }
    FunctionNode& operator=(const FunctionNode<D> &n) { NOT_IMPLEMENTED_ABORT; }
    virtual ~FunctionNode() { }

    double evalf(const double *r);

    double integrate() const;
    double dotScaling(const FunctionNode<D> &ket) const;
    double dotWavelet(const FunctionNode<D> &ket) const;

    double evalScaling(const double *r) const;
    double integrateLegendre() const;
    double integrateInterpolating() const;

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
        ar & boost::serialization::base_object<MWNode<D> >(*this);
    }
};

#endif // FUNCTIONNODE_H
