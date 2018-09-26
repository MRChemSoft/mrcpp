/**
*  Basic class for representing functions in a multiwavelet
*  representation.
*/

#pragma once

#include "MWTree.h"
#include <array>

namespace mrcpp {

template<int D>
class FunctionTree: public MWTree<D> {
public:
    FunctionTree(const MultiResolutionAnalysis<D> &mra, SharedMemory *sh_mem = nullptr);
    virtual ~FunctionTree();

    void clear();

    double integrate() const;
    virtual double evalf(const double *r);
    virtual double evalf(const std::array<double, D> &r);

    void getEndValues(Eigen::VectorXd &data);
    void setEndValues(Eigen::VectorXd &data);

    void saveTree(const std::string &file);
    void loadTree(const std::string &file);

    // In place operations
    void square();
    void power(double p);
    void rescale(double c);
    void normalize();
    void add(double c, FunctionTree<D> &inp);
    void multiply(double c, FunctionTree<D> &inp);

    int getNChunks();
    int getNChunksUsed();

    FunctionNode<D> &getEndFuncNode(int i) { return static_cast<FunctionNode<D> &>(this->getEndMWNode(i)); }
    FunctionNode<D> &getRootFuncNode(int i) { return static_cast<FunctionNode<D> &>(this->rootBox.getNode(i)); }

    SerialFunctionTree<D>* getSerialFunctionTree() { return static_cast<SerialFunctionTree<D> *>(this->serialTree_p); }
    void printSerialIndices();

    const FunctionNode<D> &getEndFuncNode(int i) const { return static_cast<const FunctionNode<D> &>(this->getEndMWNode(i)); }
    const FunctionNode<D> &getRootFuncNode(int i) const { return static_cast<const FunctionNode<D> &>(this->rootBox.getNode(i)); }

protected:
    std::ostream& print(std::ostream &o);
};

}
