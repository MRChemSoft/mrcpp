#ifndef MULTIPLICATIONVECTOR_H
#define MULTIPLICATIONVECTOR_H

#include <vector>

#include "FunctionTree.h"

template<int D> class MWMultiplier;

template<int D>
class MultiplicationVector {
public:
    void push_back(FunctionTree<D> *f) {
        this->funcs.push_back(f);
    }
    void clear() {
        this->funcs.clear();
    }
    int size() const {
        return this->funcs.size();
    }
    FunctionTree<D> &getFunc(int i) {
        if (i < 0 or i >= this->funcs.size()) MSG_ERROR("Out of bounds");
        return *this->funcs[i];
    }
    friend class MWMultiplier<D>;
protected:
    std::vector<FunctionTree<D> *> funcs;

    void clean() {
        for (int n = 0; n < this->funcs.size(); n++) {
            getFunc(n).deleteGenerated();
        }
    }
};

#endif // MULTIPLICATIONVECTOR_H
