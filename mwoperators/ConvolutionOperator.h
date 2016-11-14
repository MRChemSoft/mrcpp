#ifndef CONVOLUTIONOPERATOR_H
#define CONVOLUTIONOPERATOR_H

#include "FunctionTreeVector.h"

template<int D>
class ConvolutionOperator {
public:
    ConvolutionOperator(const MultiResolutionAnalysis<D> &mra, double pr);
    virtual ~ConvolutionOperator();

    int size() const { return this->oper_exp.size(); }
    void push_back(OperatorTree *oper) { this->oper_exp.push_back(oper); } 
    void clear(bool dealloc = false) {
        if (dealloc) {
            for (int i = 0; i < oper_exp.size(); i++) {
                if (this->oper_exp[i] != 0) delete this->oper_exp[i];
            }
        }
        this->oper_exp.clear();
    }

    int getMaxBandWidth(int depth = -1) const;
    const Eigen::VectorXi &getMaxBandWidths() const { return this->bandMax; }

    void calcBandWidths(double prec);
    void clearBandWidths();

    OperatorTree &getComponent(int i);
    const OperatorTree &getComponent(int i) const;

    OperatorTree *operator[](int i) { return this->oper_exp[i]; }
    const OperatorTree *operator[](int i) const { return this->oper_exp[i]; }

protected:
    double prec;
    Eigen::VectorXi bandMax;
    MultiResolutionAnalysis<2> oper_mra;
    MultiResolutionAnalysis<1> kern_mra;
    OperatorTreeVector oper_exp;
    FunctionTreeVector<1> kern_exp;

    void initializeOperator(GreensKernel &greens_kernel);
    void clearKernel();
    void clearOperator();

    double calcMinDistance(const MultiResolutionAnalysis<D> &MRA, double epsilon) const;
    double calcMaxDistance(const MultiResolutionAnalysis<D> &MRA) const;
};

#endif // CONVOLUTIONOPERATOR_H
