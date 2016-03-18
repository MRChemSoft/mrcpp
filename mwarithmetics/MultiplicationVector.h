#ifndef MULTIPLICATIONVECTOR_H
#define MULTIPLICATIONVECTOR_H

#include <vector>

#include "FunctionTree.h"

template<int D> class MWMultiplier;

template<int D>
class MultiplicationVector {
public:
    MultiplicationVector() : coef(1.0) { }
    ~MultiplicationVector() { }
    MultiplicationVector& operator=(const MultiplicationVector<D> &vec) {
        this->clear();
        this->coef = vec.getCoef();
        for (int i = 0; i < vec.size(); i++) {
            this->push_back(vec.getFunc(i));
        }
    }

    int size() const {
        return this->funcs.size();
    }
    void setCoef(double c) {
        this->coef = c;
    }
    void push_back(FunctionTree<D> &f, double c = 1.0) {
        this->coef *= c;
        this->funcs.push_back(&f);
    }
    void clear() {
        this->coef = 1.0;
        this->funcs.clear();
    }
    double getCoef() const {
        return this->coef;
    }
    FunctionTree<D> &getFunc(int i) {
        if (i < 0 or i >= this->funcs.size()) MSG_ERROR("Out of bounds");
        return *this->funcs[i];
    }
    friend class MWMultiplier<D>;
protected:
    double coef;
    std::vector<FunctionTree<D> *> funcs;

    void clean() {
        for (int n = 0; n < this->funcs.size(); n++) {
            getFunc(n).deleteGenerated();
        }
    }
};

#endif // MULTIPLICATIONVECTOR_H
