/**
*  Basic class for representing functions in a multiwavelet
*  representation.
*/

#pragma once

#include "MWTree.h"

namespace mrcpp {

template<int D>
class FunctionTree: public MWTree<D> {
public:
    FunctionTree(const MultiResolutionAnalysis<D> &mra, SharedMemory *sh_mem = 0);
    virtual ~FunctionTree();

    void clear();

    double integrate() const;
    virtual double dot(const FunctionTree<D> &ket);
    virtual double evalf(const double *r);

    void getEndValues(Eigen::VectorXd &data);
    void setEndValues(Eigen::VectorXd &data);

    bool saveTree(const std::string &file);
    bool loadTree(const std::string &file);

    // In place operations
    void square();
    void power(double d);
    void normalize();
    void map(const RepresentableFunction<1> &func);

    FunctionTree<D>& operator *=(double c);
    FunctionTree<D>& operator *=(const FunctionTree<D> &tree);
    FunctionTree<D>& operator +=(const FunctionTree<D> &tree);
    FunctionTree<D>& operator -=(const FunctionTree<D> &tree);

    FunctionNode<D> &getEndFuncNode(int i) { return static_cast<FunctionNode<D> &>(this->getEndMWNode(i)); }
    FunctionNode<D> &getRootFuncNode(int i) { return static_cast<FunctionNode<D> &>(this->rootBox.getNode(i)); }

    SerialFunctionTree<D>* getSerialFunctionTree() { return static_cast<SerialFunctionTree<D> *>(this->serialTree_p); }

    const FunctionNode<D> &getEndFuncNode(int i) const { return static_cast<const FunctionNode<D> &>(this->getEndMWNode(i)); }
    const FunctionNode<D> &getRootFuncNode(int i) const { return static_cast<const FunctionNode<D> &>(this->rootBox.getNode(i)); }

protected:
    std::ostream& print(std::ostream &o);
};

}
