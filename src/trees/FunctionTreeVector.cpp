#include "trees/FunctionTreeVector.h"
#include "mwutils/Printer.h"

namespace mrcpp {

template<int D>
FunctionTreeVector<D>::FunctionTreeVector(const FunctionTreeVector<D> &vec) {
    this->clear();
    for (int i = 0; i < vec.size(); i++) {
        this->coefs.push_back(vec.coefs[i]);
        this->funcs.push_back(vec.funcs[i]);
    }
}

template<int D>
FunctionTreeVector<D>& FunctionTreeVector<D>::operator=(const FunctionTreeVector<D> &vec) {
    this->clear();
    for (int i = 0; i < vec.size(); i++) {
        this->coefs.push_back(vec.coefs[i]);
        this->funcs.push_back(vec.funcs[i]);
    }
    return *this;
}

template<int D>
void FunctionTreeVector<D>::push_back(double c, FunctionTree<D> *f) {
    this->coefs.push_back(c);
    this->funcs.push_back(f);
}

template<int D>
void FunctionTreeVector<D>::push_back(FunctionTree<D> *f) {
    this->coefs.push_back(1.0);
    this->funcs.push_back(f);
}

template<int D>
void FunctionTreeVector<D>::clear(bool dealloc) {
    if (dealloc) {
        for (int i = 0; i < this->funcs.size(); i++) {
            if (this->funcs[i] != 0) delete this->funcs[i];
        }
    }
    this->coefs.clear();
    this->funcs.clear();
}

template<int D>
double FunctionTreeVector<D>::getCoef(int i) const {
    if (i < 0 or i >= this->coefs.size()) MSG_ERROR("Out of bounds");
    return this->coefs[i];
}

template<int D>
FunctionTree<D> &FunctionTreeVector<D>::getFunc(int i) {
    if (this->funcs[i] == 0) MSG_ERROR("Invalid function");
    if (i < 0 or i >= this->funcs.size()) MSG_ERROR("Out of bounds");
    return *this->funcs[i];
}

template<int D>
const FunctionTree<D> &FunctionTreeVector<D>::getFunc(int i) const {
    if (this->funcs[i] == 0) MSG_ERROR("Invalid function");
    if (i < 0 or i >= this->funcs.size()) MSG_ERROR("Out of bounds");
    return *this->funcs[i];
}

template<int D>
FunctionTree<D> *FunctionTreeVector<D>::operator[](int i) {
    if (i < 0 or i >= this->funcs.size()) MSG_ERROR("Out of bounds");
    return this->funcs[i];
}

template<int D>
const FunctionTree<D> *FunctionTreeVector<D>::operator[](int i) const {
    if (i < 0 or i >= this->funcs.size()) MSG_ERROR("Out of bounds");
    return this->funcs[i];
}

template class FunctionTreeVector<1>;
template class FunctionTreeVector<2>;
template class FunctionTreeVector<3>;

} //namespace mrcpp
