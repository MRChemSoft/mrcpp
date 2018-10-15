
/*
 *
 */

#pragma once

#include "FunctionNode.h"

namespace mrcpp {

template<int D>
class GenNode final : public FunctionNode<D> {
public:
    double getWaveletNorm() const { return 0.0; }

    void createChildren();
    void genChildren();
    void cvTransform(int kind);
    void mwTransform(int kind);

    void setValues(const Eigen::VectorXd &vec);
    void getValues(Eigen::VectorXd &vec);

    friend class SerialFunctionTree<D>;

protected:
    GenNode() : FunctionNode<D>() { }
    GenNode(const GenNode<D> &node) = delete;
    GenNode<D> &operator=(const GenNode<D> &node) = delete;
    ~GenNode() { assert(this->tree == 0); }

    double calcComponentNorm(int i) const;
    void dealloc();
    void reCompress();
};

}
