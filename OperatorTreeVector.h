#ifndef OPERATORTREEVECTOR_H
#define OPERATORTREEVECTOR_H

#include <vector>

#include "OperatorTree.h"

class OperatorTreeVector {
public:
    OperatorTreeVector() { }
    virtual ~OperatorTreeVector() { }

    int size() const {
        return this->operComp.size();
    }
    void push_back(OperatorTree &o) {
        this->operComp.push_back(&o);
    }
    void clear() {
        this->operComp.clear();
    }
    OperatorTree &getComponent(int i) {
        if (i < 0 or i >= this->operComp.size()) MSG_ERROR("Out of bounds");
        return *this->operComp[i];
    }
    const OperatorTree &getComponent(int i) const {
        if (i < 0 or i >= this->operComp.size()) MSG_ERROR("Out of bounds");
        return *this->operComp[i];
    }
protected:
    std::vector<OperatorTree *> operComp;
};


#endif // OPERATORTREEVECTOR_H
