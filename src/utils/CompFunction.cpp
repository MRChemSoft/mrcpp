#include <fstream>
#include "Printer.h"
#include "parallel.h"
#include "treebuilders/project.h"
#include "treebuilders/add.h"
#include "treebuilders/multiply.h"
#include "CompFunction.h"
#include "ComplexFunction.h"

namespace mrcpp {

  template <int D>
  MultiResolutionAnalysis<D> *defaultCompMRA = nullptr; // Global MRA

  template <int D, typename T>
  template <int D_, typename std::enable_if<D_ == 3, int>::type>
  CompFunction<D, T>::CompFunction(ComplexFunction cplxfunc){
      Ncomp = 1;
      if (std::is_same<T, ComplexDouble>::value) {
          isreal = 0;
          iscomplex = 1;
      } else {
          isreal = 1;
          iscomplex = 0;
      }
      defaultCompMRA<3> = cplxfunc.funcMRA;
      //we always copy real part
      Comp[0] = new FunctionTree<D, T>(*cplxfunc.funcMRA);
      if (not cplxfunc.hasReal()) MSG_ABORT("Input function has no real part");
      FunctionTree<D, T>::deep_copy(Comp[0], cplxfunc.real());
      if ( iscomplex ){
          //We add the imaginary part, if it exist in input function
          if (cplxfunc.hasImag()){
              ComplexDouble c;
              if(cplxfunc.conjugate())  MSG_ERROR("onjugaison not implemented");
              Comp[0].add_inplace(1.0, cplxfunc.imag());
          }
      } else if (cplxfunc.hasImag()) MSG_WARN("Complex part is truncated")
      // set metadata
      data.n1[0] = cplxfunc.spin();
      data.n2[0] = cplxfunc.occ();

      rank = cplxfunc.getRank();
  }

  template <int D, typename T>
  CompFunction<D, T>::CompFunction()
  { if (std::is_same<T, ComplexDouble>::value) {
          isreal = 0;
          iscomplex = 1;
      } else {
          isreal = 1;
          iscomplex = 0;
      }
    data.Ncomp = 0;
    Comp[0]=nullptr;
    Comp[1]=nullptr;
    Comp[2]=nullptr;
    Comp[3]=nullptr;
  }

  template <int D, typename T>
  CompFunction<D, T>::CompFunction(MultiResolutionAnalysis<D> &mra)
  { if (std::is_same<T, ComplexDouble>::value) {
          isreal = 0;
          iscomplex = 1;
      } else {
          isreal = 1;
          iscomplex = 0;
      }
    defaultCompMRA<D> = &mra;
    data.Ncomp = 0;
    Comp[0]=nullptr;
    Comp[1]=nullptr;
    Comp[2]=nullptr;
    Comp[3]=nullptr;
  }

/** @brief Copy constructor
 *
 * Shallow copy: meta data is copied along with the component pointers,
 * NO transfer of ownership.
 */
  template <int D, typename T>
  CompFunction<D, T>::CompFunction(const CompFunction<D, T> &compfunc) {
      data = compfunc.data;
      Comp[0] = compfunc.Comp[0];
      Comp[1] = compfunc.Comp[1];
      Comp[2] = compfunc.Comp[2];
      Comp[3] = compfunc.Comp[3];
  }

/** @brief Copy constructor
 *
 * Shallow copy: meta data is copied along with the component pointers,
 * NO transfer of ownership.
 */
  template <int D, typename T>
  CompFunction<D, T>::CompFunction(CompFunction<D, T> && compfunc) {
      data = compfunc.data;
      Comp[0] = compfunc.Comp[0];
      Comp[1] = compfunc.Comp[1];
      Comp[2] = compfunc.Comp[2];
      Comp[3] = compfunc.Comp[3];
  }

  template <int D, typename T>
  CompFunction<D, T> &CompFunction<D, T>::operator=(const CompFunction<D, T> &func) {
      if (this != &func) {
          this->data = func.data;
          for (int i = 0; i < Ncomp; i++) {
              this->Comp[i] = func.Comp[i];
          }
      }
      return *this;
  }

    template <int D, typename T>
    template <int D_, typename std::enable_if<D_ == 3, int>::type>
    CompFunction<D, T>::operator ComplexFunction() const {
        return ComplexFunction(*this); // const conversion
    }
    //    template <int D, typename T>
   //    template <int D_, typename std::enable_if<D_ == 3, int>::type>
   //    CompFunction<D, T>::operator ComplexFunction() {
   //        return ComplexFunction(std::move(*this)); // non-const conversion
   //    }
  //

    template CompFunction<3,double>::operator ComplexFunction() const;
//template CompFunction<3,double>::operator ComplexFunction();

    template <int D, typename T>
    void CompFunction<D, T>::flushFuncData() {
      for (int i = 0; i < Ncomp; i++) {
          Nchunks[i] = Comp[i]->getNChunksUsed();
      }
      for (int i = Ncomp; i < 4; i++) Nchunks[i] = 0;
    }
  template <int D, typename T>
  double CompFunction<D, T>::norm() const {
     double norm = squaredNorm();
     for (int i = 0; i < Ncomp; i++) {
          norm += Comp[i]->getSquareNorm();
     }
     if (norm > 0.0) norm = std::sqrt(norm);
     return norm;
  }
  template <int D, typename T>
  double CompFunction<D, T>::squaredNorm() const {
     double norm = squaredNorm();
     for (int i = 0; i < Ncomp; i++) {
          norm += Comp[i]->getSquareNorm();
     }
     return norm;
  }
  template <int D, typename T>
  void CompFunction<D, T>::alloc(int i) {
      if (defaultCompMRA<D> == nullptr) MSG_ABORT("Default MRA not yet defined");
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

/** @brief In place multiply with scalar. Fully in-place.*/
template <int D, typename T>
void CompFunction<D, T>::rescale(T c) {
    bool need_to_rescale = not(isShared()) or mpi::share_master();
    if (need_to_rescale) {
        for (int i = 0; i < Ncomp; i++) {
            Comp[i]->rescale(c);
        }
    } else MSG_ERROR("Not implemented");
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


namespace compfunc {


/** @brief Deep copy
 *
 * Deep copy: meta data is copied along with the content of each component.
 */
  template <int D, typename T>
  void deep_copy(CompFunction<D, T> *out, const CompFunction<D, T> &inp) {
      out->data = inp.data;
      for (int i = 0; i < inp.Ncomp; i++) {
          delete out->Comp[i];
          inp.Comp[i]->deep_copy(out->Comp[i]);
      }
  }


/** @brief out = a*inp_a + b*inp_b
 *
 * Recast into linear_combination.
 *
 */
template <int D, typename T>
void add(CompFunction<D, T> &out, T a, CompFunction<D, T> inp_a, T b, CompFunction<D, T> inp_b, double prec) {
    std::vector<T> coefs(2);
    coefs[0] = a;
    coefs[1] = b;

    std::vector<CompFunction<D, T>> funcs; // NB: not a CompFunctionVector, because not run in parallel!
    funcs.push_back(inp_a);
    funcs.push_back(inp_b);

    linear_combination(out, coefs, funcs, prec);
}

/** @brief out = c_0*inp_0 + c_1*inp_1 + ... + c_N*inp_N
 *
 * OMP parallel, but not MPI parallel
 */
template <int D, typename T>
    void linear_combination(CompFunction<D, T> &out, const std::vector<T> &c, std::vector<CompFunction<D, T>> &inp, double prec) {
    double thrs = MachineZero;
    bool need_to_add = not(out.isShared()) or mpi::share_master();
    for (int comp = 0; comp < inp[0].Ncomp; comp++) {
        FunctionTreeVector<D, T> fvec; // one component vector
        for (int i = 0; i < inp.size(); i++) {
            if (std::norm(c[i]) < thrs) continue;
            if (out.iscomplex and inp[i].data.conj) MSG_ERROR("conjugaison not implemented");
            fvec.push_back(std::make_tuple(c[i], inp[i].Comp[comp]));
        }
        if (need_to_add) {
            if (fvec.size() > 0) {
                if (prec < 0.0) {
                    build_grid(out.real(), fvec);
                    mrcpp::add(prec, *out.Comp[comp], fvec, 0);
                } else {
                    mrcpp::add(prec, *out.Comp[comp], fvec);
                }
            } else if (out.hasReal()) {
                out.Comp[comp]->setZero();
            }
        }
        mpi::share_function(out, 0, 9911, mpi::comm_share);
    }
}

/** @brief out = inp_a * inp_b
 *
 */
template <int D, typename T>
void multiply(CompFunction<D, T> &out, CompFunction<D, T> inp_a, CompFunction<D, T> inp_b, double prec, bool absPrec, bool useMaxNorms) {
    bool need_to_multiply = not(out.isShared()) or mpi::share_master();
    for (int comp = 0; comp < inp_a[0].Ncomp; comp++) {
        delete out.Comp[comp];
        FunctionTree<3, T> *tree = new FunctionTree<3, T>(inp_a.Comp[0].getMRA());
        T coef = 1.0;
         if (need_to_multiply) {
             if (out.iscomplex and inp_a.data.conj) MSG_ERROR("conjugaison not implemented");
             if (out.iscomplex and inp_b.data.conj) MSG_ERROR("conjugaison not implemented");
             if (prec < 0.0) {
                 // Union grid
                 build_grid(*tree, inp_a.Comp[comp]);
                 build_grid(*tree, inp_b.Comp[comp]);
                 mrcpp::multiply(prec, *tree, coef, *inp_a.Comp[comp], *inp_b.Comp[comp], 0);
             } else {
                // Adaptive grid
                 mrcpp::multiply(prec, *tree, coef, *inp_a.Comp[comp], *inp_b.Comp[comp], -1, absPrec, useMaxNorms);
             }
         }
         out.Comp[comp] = tree;
    }
    mpi::share_function(out, 0, 9911, mpi::comm_share);

}

/** @brief out = inp_a * f
 *
 *  each component is multiplied
 */
template <int D, typename T>
void multiply(CompFunction<D, T> &out, CompFunction<D, T> &inp_a, RepresentableFunction<D, T> &f, double prec, int nrefine) {
    MSG_ERROR("Not implemented");
}

/** @brief out = inp_a * f
 *
 */
template <int D, typename T>
void multiply(CompFunction<D, T>, FunctionTree<D, T> &inp_a, RepresentableFunction<D, T> &f, double prec, int nrefine) {
    MSG_ERROR("Not implemented");
}


/** @brief Compute <bra|ket> = int bra^\dag(r) * ket(r) dr.
 *
 *  Sum of component dots.
 *  Notice that the <bra| position is already complex conjugated.
 *
 */
template <int D, typename T>
T compfunc::dot(CompFunction<D, T> bra, CompFunction<D, T> ket) {
    T dotprod = 0.0;
    if (bra.data.conj or ket.data.conj) MSG_ERROR("dot with conjugaison not implemented");
    for (int comp = 0; comp < bra.Ncomp; comp++) {
        dotprod += mrcpp::dot(bra.Comp[comp], ket.Comp[comp]);
    }
    return dotprod;
}


} // namespace compfunc

} // namespace mrcpp
