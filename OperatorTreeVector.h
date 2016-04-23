#ifndef OPERATORTREEVECTOR_H
#define OPERATORTREEVECTOR_H

#include <vector>

#include "OperatorTree.h"

class OperatorTreeVector {
public:
    OperatorTreeVector() { }
    virtual ~OperatorTreeVector() { }

    int size() const { return this->operComp.size(); }

    int getMaxBandWidth(int depth = -1) const;
    void calcBandWidths(double prec);

    void push_back(OperatorTree &o) { this->operComp.push_back(&o); }
    void clear() { this->operComp.clear(); }

    OperatorTree &getComponent(int i) { return *this->operComp[i]; }
    const OperatorTree &getComponent(int i) const { return *this->operComp[i]; }

    OperatorTree *operator[](int i) { return this->operComp[i]; }
    const OperatorTree *operator[](int i) const { return this->operComp[i]; }

protected:
    Eigen::VectorXi bandMax;
    std::vector<OperatorTree *> operComp;

};


#endif // OPERATORTREEVECTOR_H
