#ifndef ADDITIONVECTOR_H
#define ADDITIONVECTOR_H

#include <vector>

#include "FunctionTree.h"

template<int D> class MWAdder;

template<int D>
class AdditionVector {
public:
    void push_back(double c, FunctionTree<D> *f) {
        this->coefs.push_back(c);
        this->funcs.push_back(f);
    }
    void clear() {
        this->coefs.clear();
        this->funcs.clear();
    }
    int size() const {
        return this->funcs.size();
    }
    double getCoef(int i) const {
        if (i < 0 or i >= this->coefs.size()) MSG_ERROR("Out of bounds");
        return this->coefs[i];
    }
    FunctionTree<D> &getFunc(int i) {
        if (i < 0 or i >= this->funcs.size()) MSG_ERROR("Out of bounds");
        return *this->funcs[i];
    }
    friend class MWAdder<D>;
protected:
    std::vector<double> coefs;
    std::vector<FunctionTree<D> *> funcs;

    void clean() {
        for (int n = 0; n < this->funcs.size(); n++) {
            getFunc(n).deleteGenerated();
        }
    }
};

#endif // ADDITIONVECTOR_H
