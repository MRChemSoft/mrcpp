#include <fstream>
#include "Printer.h"
#include "CompFunction.h"
#include "ComplexFunction.h"

namespace mrcpp {

  template <int D>
  MultiResolutionAnalysis<D> *defaultCompMRA; // Global MRA

  template <int D, typename T>
  template <int D_, typename std::enable_if<D_ == 3, int>::type>
  CompFunction<D, T>::CompFunction(T value, ComplexFunction cplxfunc)
      : Ncomp(1){
      defaultCompMRA<3> = cplxfunc.funcMRA;
      //we always copy real part
      Comp[0] = new FunctionTree<D, T>(*cplxfunc.funcMRA);
      if (not cplxfunc.hasReal()) MSG_ABORT("Input funcion has nor real part");
      deep_copy(Comp[0], cplxfunc.real());
      if (std::is_same<T, ComplexDouble>::value){
          //We add the imaginary part, if it exist in input function
          if (cplxfunc.hasImag()){
              ComplexDouble c;
              if(cplxfunc.conjugate()) c = {0.0, -1.0};
              else c = {0.0, 1.0};
              Comp[0].add_inplace(c, cplxfunc.imag());
          }
      }
      // set metadata
      data.n1[0] = cplxfunc.spin();
      data.n2[0] = cplxfunc.occ();

      rank = cplxfunc.getRank();
  }
  template <int D, typename T>
  CompFunction<D, T>::CompFunction(MultiResolutionAnalysis<D> &mra)
  { defaultCompMRA<D> = &mra;
    data.Ncomp = 0;
    Comp[0]=nullptr;
    Comp[1]=nullptr;
    Comp[2]=nullptr;
    Comp[3]=nullptr;
  }
  template <int D, typename T>
  double CompFunction<D, T>::norm() {
     double norm = squaredNorm();
     for (int i = 0; i < Ncomp; i++) {
          norm += Comp[i]->getSquareNorm();
     }
     if (norm > 0.0) norm = std::sqrt(norm);
     return norm;
  }
  template <int D, typename T>
  double CompFunction<D, T>::squaredNorm() {
     double norm = squaredNorm();
     for (int i = 0; i < Ncomp; i++) {
          norm += Comp[i]->getSquareNorm();
     }
     return norm;
  }
  template <int D, typename T>
  void CompFunction<D, T>::alloc(int i) {
      Comp[i] = new FunctionTree<D, T> (*defaultCompMRA<D>);
  }

/** @brief In place addition.
 *
 * Output is extended to union grid.
 *
 */
template <int D, typename T>
void CompFunction<D, T>::add(T c, CompFunction<D, T> inp) {
    for (int i = 0; i < Ncomp; i++) {
        if (i >= inp.Ncomp) break;
        Comp[i]->add_inplace(c,*inp.Comp[i]);
    }
    for (int i = Ncomp; i < inp.Ncomp; i++) {
        alloc(i);
        Comp[i]->add_inplace(c,*inp.Comp[i]);
    }
}

template <int D, typename T>
int CompFunction<D, T>::crop(double prec) {
    if (prec < 0.0) return 0;
    int nChunksremoved = 0;
    for (int i = 0; i < Ncomp; i++) {
        nChunksremoved += Comp[i]->crop(prec, 1.0, false);
    }
    return nChunksremoved;
}


template class  MultiResolutionAnalysis<1>;
template class  MultiResolutionAnalysis<2>;
template class  MultiResolutionAnalysis<3>;
template class CompFunction<1, double>;
template class CompFunction<2, double>;
template class CompFunction<3, double>;

template class CompFunction<1, ComplexDouble>;
template class CompFunction<2, ComplexDouble>;
template class CompFunction<3, ComplexDouble>;

} // namespace mrcpp
