#pragma once

#include <vector>

#include "FunctionTree.h"

template<int D>
class FunctionTreeVector {
public:
    FunctionTreeVector() { }
    virtual ~FunctionTreeVector() { }

    FunctionTreeVector(const FunctionTreeVector<D> &vec);
    FunctionTreeVector& operator=(const FunctionTreeVector<D> &vec);

    int size() const { return this->funcs.size(); }
    void clear(bool dealloc = false);

    void push_back(double c, FunctionTree<D> *f);
    void push_back(FunctionTree<D> *f);

    double getCoef(int i) const;
    FunctionTree<D> &getFunc(int i);
    const FunctionTree<D> &getFunc(int i) const;

    FunctionTree<D> *operator[](int i);
    const FunctionTree<D> *operator[](int i) const;

protected:
    std::vector<double> coefs;
    std::vector<FunctionTree<D> *> funcs;
};

