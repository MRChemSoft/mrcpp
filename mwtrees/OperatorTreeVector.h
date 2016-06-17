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
    const Eigen::VectorXi &getMaxBandWidths() const { return this->bandMax; }

    void calcBandWidths(double prec);
    void clearBandWidths();

    void push_back(OperatorTree *o) { this->operComp.push_back(o); }
    void clear(bool dealloc = false) {
        if (dealloc) {
            for (int i = 0; i < this->operComp.size(); i++) {
                if (this->operComp[i] != 0) delete this->operComp[i];
            }
        }
        this->operComp.clear();
    }

    OperatorTree &getComponent(int i) {
        if (this->operComp[i] == 0) MSG_ERROR("Invalid component");
        if (i < 0 or i >= this->operComp.size()) MSG_ERROR("Out of bounds");
        return *this->operComp[i];
    }
    const OperatorTree &getComponent(int i) const {
        if (this->operComp[i] == 0) MSG_ERROR("Invalid component");
        if (i < 0 or i >= this->operComp.size()) MSG_ERROR("Out of bounds");
        return *this->operComp[i];
    }

    OperatorTree *operator[](int i) { return this->operComp[i]; }
    const OperatorTree *operator[](int i) const { return this->operComp[i]; }

protected:
    Eigen::VectorXi bandMax;
    std::vector<OperatorTree *> operComp;

};


#endif // OPERATORTREEVECTOR_H
