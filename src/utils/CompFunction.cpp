#include <fstream>
#include "Printer.h"
#include "parallel.h"
#include "Bank.h"
#include "treebuilders/grid.h"
#include "trees/FunctionNode.h"
#include "treebuilders/project.h"
#include "treebuilders/add.h"
#include "treebuilders/multiply.h"
#include "CompFunction.h"

namespace mrcpp {

  template <int D>
  MultiResolutionAnalysis<D> *defaultCompMRA = nullptr; // Global MRA

  template <int D>
  CompFunction<D>::CompFunction(MultiResolutionAnalysis<D> &mra)
  { defaultCompMRA<D> = &mra;
    data.Ncomp = 0;
    for (int i = 0; i < 4; i++) CompD[i] = nullptr;
    for (int i = 0; i < 4; i++) CompC[i] = nullptr;
  }

  template <int D>
  CompFunction<D>::CompFunction()
  { data.Ncomp = 0;
    for (int i = 0; i < 4; i++) CompD[i] = nullptr;
    for (int i = 0; i < 4; i++) CompC[i] = nullptr;
  }

/*
 * Empty functions (no components defined)
 */
  template <int D>
  CompFunction<D>::CompFunction(int n1)
  { data.Ncomp = 0;
      data.n1[0] = n1;
      data.n2[0] = -1;
      data.n3[0] = 0;
      rank = 0;
      isreal = 1;
      iscomplex = 0;
      data.shared = false;
      //      if (defaultCompMRA<D> == nullptr) MSG_ABORT("Default MRA not yet defined");
      //CompD[i] = new FunctionTree<D, double> (*defaultCompMRA<D>);
    for (int i = 0; i < 4; i++) CompD[i] = nullptr;
    for (int i = 0; i < 4; i++) CompC[i] = nullptr;

  }

/*
 * Empty functions (no components defined)
 */
  template <int D>
  CompFunction<D>::CompFunction(int n1, bool share)
  { data.Ncomp = 0;
      data.n1[0] = n1;
      data.n2[0] = -1;
      data.n3[0] = 0;
      rank = 0;
      isreal = 1;
      iscomplex = 0;
      data.shared = share;
     if (share) MSG_ABORT("Not yet implemented");
      //CompD[i] = new FunctionTree<D, double> (*defaultCompMRA<D>);
    for (int i = 0; i < 4; i++) CompD[i] = nullptr;
    for (int i = 0; i < 4; i++) CompC[i] = nullptr;

  }

/** @brief Copy constructor
 *
 * Shallow copy: meta data is copied along with the component pointers,
 * NO transfer of ownership.
 */
  template <int D>
  CompFunction<D>::CompFunction(const CompFunction<D> &compfunc) {
      data = compfunc.data;
    for (int i = 0; i < 4; i++) CompD[i] = compfunc.CompD[i];
    for (int i = 0; i < 4; i++) CompC[i] = compfunc.CompC[i];
  }

/** @brief Copy constructor
 *
 * Shallow copy: meta data is copied along with the component pointers,
 * NO transfer of ownership.
 */
  template <int D>
  CompFunction<D>::CompFunction(CompFunction<D> && compfunc) {
      data = compfunc.data;
      for (int i = 0; i < 4; i++) CompD[i]=compfunc.CompD[i];
      for (int i = 0; i < 4; i++) CompC[i]=compfunc.CompC[i];
  }

  template <int D>
  CompFunction<D> &CompFunction<D>::operator=(const CompFunction<D> &compfunc) {
      if (this != &compfunc) {
          this->data = compfunc.data;
          for (int i = 0; i < Ncomp; i++) {
              CompD[i] = compfunc.CompD[i];
              CompC[i] = compfunc.CompC[i];
          }
      }
      return *this;
  }

//    template <int D>
//    template <int D_, typename std::enable_if<D_ == 3, int>::type>
//    CompFunction<D>::operator ComplexFunction() const {
 //       return ComplexFunction(*this); // const conversion
 //   }
    //    template <int D, typename T>
   //    template <int D_, typename std::enable_if<D_ == 3, int>::type>
   //    CompFunction<D, T>::operator ComplexFunction() {
   //        return ComplexFunction(std::move(*this)); // non-const conversion
   //    }
  //

 //   template CompFunction<3>::operator ComplexFunction() const;

    template <int D>
    void CompFunction<D>::flushFuncData() {
      for (int i = 0; i < Ncomp; i++) {
          if (isreal) {
              Nchunks[i] = CompD[i]->getNChunksUsed();
          } else {
              Nchunks[i] = CompC[i]->getNChunksUsed();
          }
      }
      for (int i = Ncomp; i < 4; i++) Nchunks[i] = 0;
    }
  template <int D>
  double CompFunction<D>::norm() const {
     double norm = squaredNorm();
     for (int i = 0; i < Ncomp; i++) {
          if (isreal) {
              norm += CompD[i]->getSquareNorm();
          } else {
              norm += CompC[i]->getSquareNorm();
          }
     }
     if (norm > 0.0) norm = std::sqrt(norm);
     return norm;
  }
  template <int D>
  double CompFunction<D>::squaredNorm() const {
     double norm = squaredNorm();
     for (int i = 0; i < Ncomp; i++) {
          if (isreal) {
              norm += CompD[i]->getSquareNorm();
          } else {
              norm += CompC[i]->getSquareNorm();
          }
     }
     return norm;
  }
  template <int D>
  void CompFunction<D>::alloc(int i) {
      if (defaultCompMRA<D> == nullptr) MSG_ABORT("Default MRA not yet defined");
      if (CompD[i] != nullptr) delete CompD[i];
      if (CompC[i] != nullptr) delete CompC[i];
      if (isreal) {
          CompD[i] = new FunctionTree<D, double> (*defaultCompMRA<D>);
      } else {
          CompC[i] = new FunctionTree<D, ComplexDouble> (*defaultCompMRA<D>);
      }
      Ncomp = std::max(Ncomp, i + 1);
 }

template <int D>
void CompFunction<D>::free() {
    //TODO: shared memory handling
    for (int i = 0; i < Ncomp; i++) {
        if (CompD[i]!= nullptr) {
            delete CompD[i];
        }
        if (CompC[i]!= nullptr) {
            delete CompC[i];
        }
    }
}

template <int D>
int CompFunction<D>::getSizeNodes() const {
    int size_mb = 0; // Memory size in kB
    for (int i = 0; i < Ncomp; i++) {
        if (CompD[i]!= nullptr) size_mb +=CompD[i]->getSizeNodes();
        if (CompC[i]!= nullptr) size_mb +=CompC[i]->getSizeNodes();
    }
    return size_mb;
}

template <int D>
int CompFunction<D>::getNNodes() const {
    int nNodes = 0;
     for (int i = 0; i < Ncomp; i++) {
        if (CompD[i]!= nullptr) nNodes +=CompD[i]->getSizeNodes();
        if (CompC[i]!= nullptr) nNodes +=CompC[i]->getSizeNodes();
    }
    return nNodes;
}

template <int D>
FunctionTree<D, double> &CompFunction<D>::real(int i) {
    if (!isreal) MSG_ABORT("not real function");
    return *CompD[i];
}
template <int D> //NB: should return CompC in the future
FunctionTree<D, double>  &CompFunction<D>::imag(int i) {
    if (!iscomplex) MSG_ABORT("not complex function");
    return *CompD[i];
}

template <int D>
const FunctionTree<D, double> &CompFunction<D>::real(int i) const {
    if (!isreal) MSG_ABORT("not real function");
    return *CompD[i];
}
template <int D> //NB: should return CompC in the future
const FunctionTree<D, double> &CompFunction<D>::imag(int i) const {
    if (!iscomplex) MSG_ABORT("not complex function");
    return *CompD[i];
}

 /* for backwards compatibility */
template <int D>
void CompFunction<D>::setReal(FunctionTree<D, double> *tree, int i) {
      if (CompD[i] != nullptr) delete CompD[i];
      if (iscomplex) MSG_ERROR("cannot write real tree into complex function");
      CompD[i] = tree;
      if (tree != nullptr) {
          isreal = 1;
          Ncomp = std::max(Ncomp, i + 1);
      } else {Ncomp = std::min(Ncomp, i);}
}
    /*
template <int D>
void CompFunction<D>::set(FunctionTree<D, ComplexDouble> *tree, int i) {
      if (CompC[i] != nullptr) delete CompD[i];
      if (isreal) MSG_ERROR("cannot write comlex tree into complex function");
      CompC[i] = tree;
      if (tree != nullptr) {
          iscomplex = 1;
          Ncomp = std::max(Ncomp, i + 1);
      } else {Ncomp = std::min(Ncomp, i);}
      } */

/** @brief In place addition.
 *
 * Output is extended to union grid.
 *
 */
template <int D>
void CompFunction<D>::add(ComplexDouble c, CompFunction<D> inp) {
    for (int i = 0; i < inp.Ncomp; i++) {
        if (i >= Ncomp) alloc(i);
        if (isreal) {
            CompD[i]->add_inplace(c.real(),*inp.CompD[i]);
        } else {
            CompC[i]->add_inplace(c,*inp.CompC[i]);
        }
    }
}


template <int D>
int CompFunction<D>::crop(double prec) {
    if (prec < 0.0) return 0;
    int nChunksremoved = 0;
    for (int i = 0; i < Ncomp; i++) {
        if (isreal) {
            nChunksremoved += CompD[i]->crop(prec, 1.0, false);
        } else {
            nChunksremoved += CompC[i]->crop(prec, 1.0, false);
        }
    }
    return nChunksremoved;
}

/** @brief In place multiply with scalar. Fully in-place.*/
template <int D>
void CompFunction<D>::rescale(ComplexDouble c) {
    bool need_to_rescale = not(isShared()) or mpi::share_master();
    if (need_to_rescale) {
        for (int i = 0; i < Ncomp; i++) {
            if (isreal) {
                CompD[i]->rescale(c.real());
            } else {
                CompC[i]->rescale(c);
            }
        }
    } else MSG_ERROR("Not implemented");
}


template class  MultiResolutionAnalysis<1>;
template class  MultiResolutionAnalysis<2>;
template class  MultiResolutionAnalysis<3>;
template class CompFunction<1>;
template class CompFunction<2>;
template class CompFunction<3>;

/** @brief Deep copy
 *
 * Deep copy: meta data is copied along with the content of each component.
 */
  template <int D>
  void deep_copy(CompFunction<D> *out, const CompFunction<D> &inp) {
      out->data = inp.data;
      for (int i = 0; i < inp.Ncomp; i++) {
          if (inp.isreal) {
              delete out->CompD[i];
              inp.CompD[i]->deep_copy(out->CompD[i]);
          } else {
              delete out->CompC[i];
              inp.CompC[i]->deep_copy(out->CompC[i]);
          }
      }
  }


/** @brief Deep copy
 *
 * Deep copy: meta data is copied along with the content of each component.
 */
  template <int D>
  void deep_copy(CompFunction<D> &out, const CompFunction<D> &inp) {
      out.data = inp.data;
      for (int i = 0; i < inp.Ncomp; i++) {
          if (inp.isreal) {
              delete out.CompD[i];
              inp.CompD[i]->deep_copy(out.CompD[i]);
          } else {
              delete out.CompC[i];
              inp.CompC[i]->deep_copy(out.CompC[i]);
          }
      }
  }

/** @brief out = a*inp_a + b*inp_b
 *
 * Recast into linear_combination.
 *
 */
template <int D>
void add(CompFunction<D> &out, ComplexDouble a, CompFunction<D> inp_a, ComplexDouble b, CompFunction<D> inp_b, double prec) {
    std::vector<ComplexDouble> coefs(2);
    coefs[0] = a;
    coefs[1] = b;

    std::vector<CompFunction<D>> funcs; // NB: not a CompFunctionVector, because not run in parallel!
    funcs.push_back(inp_a);
    funcs.push_back(inp_b);

    linear_combination(out, coefs, funcs, prec);
}

/** @brief out = c_0*inp_0 + c_1*inp_1 + ... + c_N*inp_N
 *
 * OMP parallel, but not MPI parallel
 */
template <int D>
    void linear_combination(CompFunction<D> &out, const std::vector<ComplexDouble> &c, std::vector<CompFunction<D>> &inp, double prec) {
    double thrs = MachineZero;
    bool need_to_add = not(out.isShared()) or mpi::share_master();
    for (int comp = 0; comp < inp[0].Ncomp; comp++) {
        if (inp[0].isreal) {
            FunctionTreeVector<D, double> fvec; // one component vector
            for (int i = 0; i < inp.size(); i++) {
                if (std::norm(c[i]) < thrs) continue;
                fvec.push_back(std::make_tuple(c[i].real(), inp[i].CompD[comp]));
            }
            if (need_to_add) {
                if (fvec.size() > 0) {
                    if (prec < 0.0) {
                        build_grid(out.CompD[comp], fvec);
                        mrcpp::add(prec, *out.CompD[comp], fvec, 0);
                    } else {
                        mrcpp::add(prec, *out.CompD[comp], fvec);
                    }
                }
            } else if (out.hasReal()) {
                out.CompD[comp]->setZero();
            }
        } else {
            FunctionTreeVector<D, ComplexDouble> fvec; // one component vector
            for (int i = 0; i < inp.size(); i++) {
                if (std::norm(c[i]) < thrs) continue;
                if (inp[i].data.conj) MSG_ERROR("conjugaison not implemented");
                fvec.push_back(std::make_tuple(c[i], inp[i].CompD[comp]));
            }
            if (need_to_add) {
                if (fvec.size() > 0) {
                    if (prec < 0.0) {
                        build_grid(out.CompC[comp], fvec);
                        mrcpp::add(prec, *out.CompC[comp], fvec, 0);
                    } else {
                        mrcpp::add(prec, *out.CompC[comp], fvec);
                    }
                }
            } else if (out.hasReal()) {
                out.CompC[comp]->setZero();
            }
        }
        mpi::share_function(out, 0, 9911, mpi::comm_share);
    }
}

/** @brief out = inp_a * inp_b
 *
 */
template <int D>
void multiply(CompFunction<D> &out, CompFunction<D> inp_a, CompFunction<D> inp_b, double prec, bool absPrec, bool useMaxNorms) {
    bool need_to_multiply = not(out.isShared()) or mpi::share_master();
    for (int comp = 0; comp < inp_a[0].Ncomp; comp++) {
        if (inp_a.isreal and inp_b.isreal) {
            delete out.CompD[comp];
            FunctionTree<D, double> *tree = new FunctionTree<D, double>(inp_a.CompD[0]->getMRA());
            double coef = 1.0;
            if (need_to_multiply) {
                if (out.iscomplex and inp_a.data.conj) MSG_ERROR("conjugaison not implemented");
                if (out.iscomplex and inp_b.data.conj) MSG_ERROR("conjugaison not implemented");
                if (prec < 0.0) {
                    // Union grid
                    build_grid(*tree, inp_a.CompD[comp]);
                    build_grid(*tree, inp_b.CompD[comp]);
                    mrcpp::multiply(prec, *tree, coef, *inp_a.CompD[comp], *inp_b.CompD[comp], 0);
                } else {
                    // Adaptive grid
                    mrcpp::multiply(prec, *tree, coef, *inp_a.CompD[comp], *inp_b.CompD[comp], -1, absPrec, useMaxNorms);
                }
            }
            out.CompD[comp] = tree;
        } else {
            delete out.CompC[comp];
            FunctionTree<D, ComplexDouble> *tree = new FunctionTree<D, ComplexDouble>(inp_a.CompC[0]->getMRA());
            ComplexDouble coef = 1.0;
            if (need_to_multiply) {
                if (out.iscomplex and inp_a.data.conj) MSG_ERROR("conjugaison not implemented");
                if (out.iscomplex and inp_b.data.conj) MSG_ERROR("conjugaison not implemented");
                if (prec < 0.0) {
                    // Union grid
                    build_grid(*tree, inp_a.CompC[comp]);
                    build_grid(*tree, inp_b.CompC[comp]);
                    mrcpp::multiply(prec, *tree, coef, *inp_a.CompC[comp], *inp_b.CompC[comp], 0);
                } else {
                    // Adaptive grid
                    mrcpp::multiply(prec, *tree, coef, *inp_a.CompC[comp], *inp_b.CompC[comp], -1, absPrec, useMaxNorms);
                }
            }
            out.CompC[comp] = tree;
        }
    }
    mpi::share_function(out, 0, 9911, mpi::comm_share);

}

/** @brief out = inp_a * f
 *
 *  each component is multiplied
 */
template <int D, typename T>
void multiply(CompFunction<D> &out, CompFunction<D> &inp_a, RepresentableFunction<D, T> &f, double prec, int nrefine) {
    MSG_ERROR("Not implemented");
}

/** @brief out = inp_a * f
 *
 */
template <int D, typename T>
void multiply(CompFunction<D>, FunctionTree<D, T> &inp_a, RepresentableFunction<D, T> &f, double prec, int nrefine) {
    MSG_ERROR("Not implemented");
}


/** @brief Compute <bra|ket> = int bra^\dag(r) * ket(r) dr.
 *
 *  Sum of component dots.
 *  Notice that the <bra| position is already complex conjugated.
 *
 */
template <int D>
ComplexDouble dot(CompFunction<D> bra, CompFunction<D> ket) {
    ComplexDouble dotprod = 0.0;
    if (bra.data.conj or ket.data.conj) MSG_ERROR("dot with conjugaison not implemented");
    for (int comp = 0; comp < bra.Ncomp; comp++) {
          if (bra.isreal and ket.isreal) {
              dotprod += mrcpp::dot(*bra.CompD[comp], *ket.CompD[comp]);
          } else {
              dotprod += mrcpp::dot(*bra.CompC[comp], *ket.CompC[comp]);
          }
    }
    if (bra.isreal and ket.isreal) {
        return dotprod.real();
    } else {
        return dotprod;
    }
}


template <int D, typename T>
void project(CompFunction<D> &out, std::function<double(const Coord<D> &r)> f, double prec) {
if (std::is_same<T, double>::value) {
    bool need_to_project = not(out.isShared()) or mpi::share_master();
    out.isreal = 1;
    out.iscomplex = 0;
    if(out.Ncomp < 1) out.alloc(0);
    if (need_to_project) mrcpp::project<D, double>(prec, out.CompD[0], f);
    mpi::share_function(out, 0, 123123, mpi::comm_share);
}
}

template <int D, typename T>
void project(CompFunction<D> &out, std::function<ComplexDouble(const Coord<D> &r)> f, double prec) {
if (std::is_same<T, ComplexDouble>::value) {
    bool need_to_project = not(out.isShared()) or mpi::share_master();
    out.isreal = 0;
    out.iscomplex = 1;
    if(out.Ncomp < 1) out.alloc(0);
    if (need_to_project) mrcpp::project<D, ComplexDouble>(prec, out.CompC[0], f);
    mpi::share_function(out, 0, 123123, mpi::comm_share);
}
 }
template <int D>
void project(CompFunction<D> &out, RepresentableFunction<D, double> &f, double prec) {
    bool need_to_project = not(out.isShared()) or mpi::share_master();
    out.isreal = 1;
    out.iscomplex = 0;
    if(out.Ncomp < 1) out.alloc(0);
    if (need_to_project) mrcpp::project<D, double>(prec, out.CompD[0], f);
    mpi::share_function(out, 0, 132231, mpi::comm_share);
}
template <int D>
void project(CompFunction<D> &out, RepresentableFunction<D, ComplexDouble> &f, double prec) {
    bool need_to_project = not(out.isShared()) or mpi::share_master();
    out.isreal = 0;
    out.iscomplex = 1;
    if(out.Ncomp < 1) out.alloc(0);
    if (need_to_project) mrcpp::project<D, ComplexDouble>(prec, out.CompC[0], f);
    mpi::share_function(out, 0, 132231, mpi::comm_share);
 }

// CompFunctionVector


CompFunctionVector::CompFunctionVector(int N):
    std::vector<CompFunction<3>*>(N) {
    for (int i = 0; i < N; i++) (*this)[i]->rank = i;
    vecMRA = defaultCompMRA<3>;
}
void CompFunctionVector::distribute() {
    for (int i = 0; i < this->size(); i++) (*this)[i]->rank = i;
}


/** @brief Make a linear combination of functions
 *
 * Uses "local" representation: treats one node at a time.
 * For each node, all functions are transformed simultaneously
 * by a dense matrix multiplication.
 * Phi input functions, Psi output functions
 *
 */
void rotate(CompFunctionVector &Phi, const ComplexMatrix &U, CompFunctionVector &Psi, double prec) {

    // The principle of this routine is that nodes are rotated one by one using matrix multiplication.
    // The routine does avoid when possible to move data, but uses pointers and indices manipulation.
    // MPI version does not use OMP yet, Serial version uses OMP
    // size of input is N, size of output is M
    int N = Phi.size();
    int M = Psi.size();
    if (U.rows() < N) MSG_ABORT("Incompatible number of rows for U matrix");
    if (U.cols() < M) MSG_ABORT("Incompatible number of columns for U matrix");

    // 1) make union tree without coefficients
    FunctionTree<3> refTree(*Phi.vecMRA);
    mpi::allreduce_Tree_noCoeff(refTree, Phi, mpi::comm_wrk);

    int sizecoeff = (1 << refTree.getDim()) * refTree.getKp1_d();
    int sizecoeffW = ((1 << refTree.getDim()) - 1) * refTree.getKp1_d();
    std::vector<double> scalefac_ref;
    std::vector<double *> coeffVec_ref; // not used!
    std::vector<int> indexVec_ref;      // serialIx of the nodes
    std::vector<int> parindexVec_ref;   // serialIx of the parent nodes
    int max_ix;
    // get a list of all nodes in union tree, identified by their serialIx indices
    refTree.makeCoeffVector(coeffVec_ref, indexVec_ref, parindexVec_ref, scalefac_ref, max_ix, refTree);
    int max_n = indexVec_ref.size();

   // 2) We work with real numbers only. Make real blocks for U matrix
    bool UhasReal = false;
    bool UhasImag = false;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            if (std::abs(U(i, j).real()) > 10*MachineZero) UhasReal = true;
            if (std::abs(U(i, j).imag()) > 10*MachineZero) UhasImag = true;
        }
    }

    IntVector PsihasReIm = IntVector::Zero(2);
    for (int j = 0; j < N; j++) {
        if (!mpi::my_func(j)) continue;
        PsihasReIm[0] = (Phi[j]->hasReal()) ? 1 : 0;
        PsihasReIm[1] = (Phi[j]->hasImag()) ? 1 : 0;
    }
    mpi::allreduce_vector(PsihasReIm, mpi::comm_wrk);
    if (not PsihasReIm[0] and not PsihasReIm[1]) {
        return; // do nothing
    }

    bool makeReal = (UhasReal and PsihasReIm[0]) or (UhasImag and PsihasReIm[1]);
    bool makeImag = (UhasReal and PsihasReIm[1]) or (UhasImag and PsihasReIm[0]);

    for (int j = 0; j < M; j++) {
        if (!mpi::my_func(j)) continue;
        if (not makeReal and Psi[j]->hasReal()) Psi[j]->free(NUMBER::Real);
        if (not makeImag and Psi[j]->hasImag()) Psi[j]->free(NUMBER::Imag);
    }

    if (not makeReal and not makeImag) { return; }

    int Neff = N;               // effective number of input orbitals
    int Meff = M;               // effective number of output orbitals
    if (makeImag) Neff = 2 * N; // Imag and Real treated independently. We always use real part of U
    if (makeImag) Meff = 2 * M; // Imag and Real treated independently. We always use real part of U

    IntVector conjMat = IntVector::Zero(Neff);
    for (int j = 0; j < Neff; j++) {
        if (!mpi::my_func(j % N)) continue;
        conjMat[j] = (Phi[j % N]->conjugate()) ? -1 : 1;
    }
    mpi::allreduce_vector(conjMat, mpi::comm_wrk);

    // we make a real matrix = U,  but organized as one or four real blocks
    // out_r = U_rr*in_r - U_ir*in_i*conjMat
    // out_i = U_ri*in_r - U_ii*in_i*conjMat
    // the first index of U is the one used on input Phi
    DoubleMatrix Ureal(Neff, Meff); // four blocks, for rr ri ir ii
    for (int j = 0; j < Neff; j++) {
        for (int i = 0; i < Meff; i++) {
            double sign = 1.0;
            if (j < N and i < M) {
                // real U applied on real Phi
                Ureal(j, i) = U.real()(j % N, i % M);
            } else if (j >= N and i >= M) {
                // real U applied on imag Phi
                Ureal(j, i) = conjMat[j] * U.real()(j % N, i % M);
            } else if (j < N and i >= M) {
                // imag U applied on real Phi
                Ureal(j, i) = U.imag()(j % N, i % M);
            } else {
                // imag U applied on imag Phi
                Ureal(j, i) = -1.0 * conjMat[j] * U.imag()(j % N, i % M);
            }
        }
    }

    // 3) In the serial case we store the coeff pointers in coeffVec. In the mpi case the coeff are stored in the bank

    bool serial = mpi::wrk_size == 1; // flag for serial/MPI switch
    BankAccount nodesPhi;             // to put the original nodes
    BankAccount nodesRotated;         // to put the rotated nodes

    // used for serial only:
    std::vector<std::vector<double *>> coeffVec(Neff);
    std::vector<std::vector<int>> indexVec(Neff);   // serialIx of the nodes
    std::map<int, std::vector<int>> node2orbVec;    // for each node index, gives a vector with the indices of the orbitals using this node
    std::vector<std::map<int, int>> orb2node(Neff); // for a given orbital and a given node, gives the node index in the
                                                    // orbital given the node index in the reference tree
    if (serial) {

        // make list of all coefficients (coeffVec), and their reference indices (indexVec)
        std::vector<int> parindexVec; // serialIx of the parent nodes
        std::vector<double> scalefac;
        for (int j = 0; j < N; j++) {
            // make vector with all coef pointers and their indices in the union grid
            if (Phi[j]->hasReal()) {
                Phi[j]->real().makeCoeffVector(coeffVec[j], indexVec[j], parindexVec, scalefac, max_ix, refTree);
                // make a map that gives j from indexVec
                int orb_node_ix = 0;
                for (int ix : indexVec[j]) {
                    orb2node[j][ix] = orb_node_ix++;
                    if (ix < 0) continue;
                    node2orbVec[ix].push_back(j);
                }
            }
            if (Phi[j]->hasImag()) {
                Phi[j]->imag().makeCoeffVector(coeffVec[j + N], indexVec[j + N], parindexVec, scalefac, max_ix, refTree);
                // make a map that gives j from indexVec
                int orb_node_ix = 0;
                for (int ix : indexVec[j + N]) {
                    orb2node[j + N][ix] = orb_node_ix++;
                    if (ix < 0) continue;
                    node2orbVec[ix].push_back(j + N);
                }
            }
        }
    } else { // MPI case

        // send own nodes to bank, identifying them through the serialIx of refTree
        save_nodes(Phi, refTree, nodesPhi);
        mpi::barrier(mpi::comm_wrk); // required for now, as the blockdata functionality has no queue yet.
    }

    // 4) rotate all the nodes
    IntMatrix split_serial;                             // in the serial case all split are stored in one array
    std::vector<std::vector<double *>> coeffpVec(Meff); // to put pointers to the rotated coefficient for each orbital in serial case
    std::vector<std::map<int, int>> ix2coef(Meff);      // to find the index in for example rotCoeffVec[] corresponding to a serialIx
    int csize;                                          // size of the current coefficients (different for roots and branches)
    std::vector<DoubleMatrix> rotatedCoeffVec;          // just to ensure that the data from rotatedCoeff is not deleted, since we point to it.
    // j indices are for unrotated orbitals, i indices are for rotated orbitals
    if (serial) {
        std::map<int, int> ix2coef_ref;   // to find the index n corresponding to a serialIx
        split_serial.resize(Meff, max_n); // not use in the MPI case
        for (int n = 0; n < max_n; n++) {
            int node_ix = indexVec_ref[n]; // SerialIx for this node in the reference tree
            ix2coef_ref[node_ix] = n;
            for (int i = 0; i < Meff; i++) split_serial(i, n) = 1;
        }

        std::vector<int> nodeReady(max_n, 0); // To indicate to OMP threads that the parent is ready (for splits)

        // assumes the nodes are ordered such that parent are treated before children. BFS or DFS ok.
        // NB: the n must be traversed approximately in right order: Thread n may have to wait until som other preceding
        // n is finished.
#pragma omp parallel for schedule(dynamic)
        for (int n = 0; n < max_n; n++) {
            int csize;
            int node_ix = indexVec_ref[n]; // SerialIx for this node in the reference tree
            // 4a) make a dense contiguous matrix with the coefficient from all the orbitals using node n
            std::vector<int> orbjVec; // to remember which orbital correspond to each orbVec.size();
            if (node2orbVec[node_ix].size() <= 0) continue;
            csize = sizecoeffW;
            if (parindexVec_ref[n] < 0) csize = sizecoeff; // for root nodes we include scaling coeff

            int shift = sizecoeff - sizecoeffW; // to copy only wavelet part
            if (parindexVec_ref[n] < 0) shift = 0;
            DoubleMatrix coeffBlock(csize, node2orbVec[node_ix].size());
            for (int j : node2orbVec[node_ix]) { // loop over indices of the orbitals using this node
                int orb_node_ix = orb2node[j][node_ix];
                for (int k = 0; k < csize; k++) coeffBlock(k, orbjVec.size()) = coeffVec[j][orb_node_ix][k + shift];
                orbjVec.push_back(j);
            }

            // 4b) make a list of rotated orbitals needed for this node
            // OMP must wait until parent is ready
            while (parindexVec_ref[n] >= 0 and nodeReady[ix2coef_ref[parindexVec_ref[n]]] == 0) {
#pragma omp flush
            };

            std::vector<int> orbiVec;
            for (int i = 0; i < Meff; i++) { // loop over all rotated orbitals
                if (not makeReal and i < M) continue;
                if (not makeImag and i >= M) continue;
                if (parindexVec_ref[n] >= 0 and split_serial(i, ix2coef_ref[parindexVec_ref[n]]) == 0) continue; // parent node has too small wavelets
                orbiVec.push_back(i);
            }

            // 4c) rotate this node
            DoubleMatrix Un(orbjVec.size(), orbiVec.size()); // chunk of U, with reorganized indices
            for (int i = 0; i < orbiVec.size(); i++) {       // loop over rotated orbitals
                for (int j = 0; j < orbjVec.size(); j++) { Un(j, i) = Ureal(orbjVec[j], orbiVec[i]); }
            }
            DoubleMatrix rotatedCoeff(csize, orbiVec.size());
            // HERE IT HAPPENS!
            rotatedCoeff.noalias() = coeffBlock * Un; // Matrix mutiplication

            // 4d) store and make rotated node pointers
            // for now we allocate in buffer, in future could be directly allocated in the final trees
            double thres = prec * prec * scalefac_ref[n] * scalefac_ref[n];
            // make all norms:
            for (int i = 0; i < orbiVec.size(); i++) {
                // check if parent must be split
                if (parindexVec_ref[n] == -1 or split_serial(orbiVec[i], ix2coef_ref[parindexVec_ref[n]])) {
                    // mark this node for this orbital for later split
#pragma omp critical
                    {
                        ix2coef[orbiVec[i]][node_ix] = coeffpVec[orbiVec[i]].size();
                        coeffpVec[orbiVec[i]].push_back(&(rotatedCoeff(0, i))); // list of coefficient pointers
                    }
                    // check norms for split
                    double wnorm = 0.0; // rotatedCoeff(k, i) is already in cache here
                    int kstart = 0;
                    if (parindexVec_ref[n] < 0) kstart = sizecoeff - sizecoeffW; // do not include scaling, even for roots
                    for (int k = kstart; k < csize; k++) wnorm += rotatedCoeff(k, i) * rotatedCoeff(k, i);
                    if (thres < wnorm or prec < 0)
                        split_serial(orbiVec[i], n) = 1;
                    else
                        split_serial(orbiVec[i], n) = 0;
                } else {
                    ix2coef[orbiVec[i]][node_ix] = max_n + 1; // should not be used
                    split_serial(orbiVec[i], n) = 0;          // do not split if parent does not need to be split
                }
            }
            nodeReady[n] = 1;
#pragma omp critical
            {
                // this ensures that rotatedCoeff is not deleted, when getting out of scope
                rotatedCoeffVec.push_back(std::move(rotatedCoeff));
            }
        }
    } else { // MPI case

        // TODO? rotate in bank, so that we do not get and put. Requires clever handling of splits.
        std::vector<double> split(Meff, -1.0);    // which orbitals need splitting (at a given node). For now double for compatibilty with bank
        std::vector<double> needsplit(Meff, 1.0); // which orbitals need splitting
        BankAccount nodeSplits;
        mpi::barrier(mpi::comm_wrk); // required for now, as the blockdata functionality has no queue yet.

        DoubleMatrix coeffBlock(sizecoeff, Neff);
        max_ix++; // largest node index + 1. to store rotated orbitals with different id
        TaskManager tasks(max_n);
        for (int nn = 0; nn < max_n; nn++) {
            int n = tasks.next_task();
            if (n < 0) break;
            double thres = prec * prec * scalefac_ref[n] * scalefac_ref[n];
            // 4a) make list of orbitals that should split the parent node, i.e. include this node
            int parentid = parindexVec_ref[n];
            if (parentid == -1) {
                // root node, split if output needed
                for (int i = 0; i < M; i++) {
                    if (makeReal)
                        split[i] = 1.0;
                    else
                        split[i] = -1.0;
                }
                for (int i = N; i < Meff; i++) {
                    if (makeImag)
                        split[i] = 1.0;
                    else
                        split[i] = -1.0;
                }
                csize = sizecoeff;
            } else {
                // note that it will wait until data is available
                nodeSplits.get_data(parentid, Meff, split.data());
                csize = sizecoeffW;
            }
            std::vector<int> orbiVec;
            std::vector<int> orbjVec;
            for (int i = 0; i < Meff; i++) {  // loop over rotated orbitals
                if (split[i] < 0.0) continue; // parent node has too small wavelets
                orbiVec.push_back(i);
            }

            // 4b) rotate this node
            DoubleMatrix coeffBlock(csize, Neff); // largest possible used size
            nodesPhi.get_nodeblock(indexVec_ref[n], coeffBlock.data(), orbjVec);
            coeffBlock.conservativeResize(Eigen::NoChange, orbjVec.size()); // keep only used part

            // chunk of U, with reorganized indices and separate blocks for real and imag:
            DoubleMatrix Un(orbjVec.size(), orbiVec.size());
            DoubleMatrix rotatedCoeff(csize, orbiVec.size());

            for (int i = 0; i < orbiVec.size(); i++) {     // loop over included rotated real and imag part of orbitals
                for (int j = 0; j < orbjVec.size(); j++) { // loop over input orbital, possibly imaginary parts
                    Un(j, i) = Ureal(orbjVec[j], orbiVec[i]);
                }
            }

            // HERE IT HAPPENS
            rotatedCoeff.noalias() = coeffBlock * Un; // Matrix mutiplication

            // 3c) find which orbitals need to further refine this node, and store rotated node (after each other while
            // in cache).
            for (int i = 0; i < orbiVec.size(); i++) { // loop over rotated orbitals
                needsplit[orbiVec[i]] = -1.0;          // default, do not split
                // check if this node/orbital needs further refinement
                double wnorm = 0.0;
                int kwstart = csize - sizecoeffW; // do not include scaling
                for (int k = kwstart; k < csize; k++) wnorm += rotatedCoeff.col(i)[k] * rotatedCoeff.col(i)[k];
                if (thres < wnorm or prec < 0) needsplit[orbiVec[i]] = 1.0;
                nodesRotated.put_nodedata(orbiVec[i], indexVec_ref[n] + max_ix, csize, rotatedCoeff.col(i).data());
            }
            nodeSplits.put_data(indexVec_ref[n], Meff, needsplit.data());
        }
        mpi::barrier(mpi::comm_wrk); // wait until all rotated nodes are ready
    }

    // 5) reconstruct trees using rotated nodes.

    // only serial case can use OMP, because MPI cannot be used by threads
    if (serial) {
        // OMP parallelized, but does not scale well, because the total memory bandwidth is a bottleneck. (the main
        // operation is writing the coefficient into the tree)

#pragma omp parallel for schedule(static)
        for (int j = 0; j < Meff; j++) {
            if (coeffpVec[j].size()==0) continue;
            if (j < M) {
                if (!Psi[j]->hasReal()) Psi[j]->alloc(0);
                Psi[j]->real().clear();
                Psi[j]->real().makeTreefromCoeff(refTree, coeffpVec[j], ix2coef[j], prec);
            } else {
                if (!Psi[j % M]->hasImag()) Psi[j % M]->alloc(0);
                Psi[j % M]->imag().clear();
                Psi[j % M]->imag().makeTreefromCoeff(refTree, coeffpVec[j], ix2coef[j], prec);
            }
        }

    } else { // MPI case

        for (int j = 0; j < Meff; j++) {
            if (not mpi::my_func(j % M)) continue;
            // traverse possible nodes, and stop descending when norm is zero (leaf in out[j])
            std::vector<double *> coeffpVec; //
            std::map<int, int> ix2coef;      // to find the index in coeffVec[] corresponding to a serialIx
            int ix = 0;
            std::vector<double *> pointerstodelete; // list of temporary arrays to clean up
            for (int ibank = 0; ibank < mpi::bank_size; ibank++) {
                std::vector<int> nodeidVec;
                double *dataVec; // will be allocated by bank
                nodesRotated.get_orbblock(j, dataVec, nodeidVec, ibank);
                if (nodeidVec.size() > 0) pointerstodelete.push_back(dataVec);
                int shift = 0;
                for (int n = 0; n < nodeidVec.size(); n++) {
                    assert(nodeidVec[n] - max_ix >= 0);                // unrotated nodes have been deleted
                    assert(ix2coef.count(nodeidVec[n] - max_ix) == 0); // each nodeid treated once
                    ix2coef[nodeidVec[n] - max_ix] = ix++;
                    csize = sizecoeffW;
                    if (parindexVec_ref[nodeidVec[n] - max_ix] < 0) csize = sizecoeff;
                    coeffpVec.push_back(&dataVec[shift]); // list of coeff pointers
                    shift += csize;
                }
            }
            if (j < M) {
                // Real part
                if (!Psi[j]->hasReal()) Psi[j]->alloc(0);
                Psi[j]->real().clear();
                Psi[j]->real().makeTreefromCoeff(refTree, coeffpVec, ix2coef, prec);
            } else {
                // Imag part
                if (!Psi[j % M]->hasImag()) Psi[j % M]->alloc(0);
                Psi[j % M]->imag().clear();
                Psi[j % M]->imag().makeTreefromCoeff(refTree, coeffpVec, ix2coef, prec);
            }
            for (double *p : pointerstodelete) delete[] p;
            pointerstodelete.clear();
        }
    }
}


void rotate(CompFunctionVector &Phi, const ComplexMatrix &U, double prec) {
    rotate(Phi, U, Phi, prec);
    return;
}

/** @brief Save all nodes in bank; identify them using serialIx from refTree
 * shift is a shift applied in the id
 */
void save_nodes(CompFunctionVector &Phi, FunctionTree<3> &refTree, BankAccount &account, int sizes) {
    int sizecoeff = (1 << refTree.getDim()) * refTree.getKp1_d();
    int sizecoeffW = ((1 << refTree.getDim()) - 1) * refTree.getKp1_d();
    int max_nNodes = refTree.getNNodes();
    std::vector<double *> coeffVec;
    std::vector<double> scalefac;
    std::vector<int> indexVec;    // SerialIx of the node in refOrb
    std::vector<int> parindexVec; // SerialIx of the parent node
    int N = Phi.size();
    int max_ix;
    for (int j = 0; j < N; j++) {
        if (not mpi::my_func(j)) continue;
        // make vector with all coef address and their index in the union grid
        if (Phi[j]->hasReal()) {
            Phi[j]->real().makeCoeffVector(coeffVec, indexVec, parindexVec, scalefac, max_ix, refTree);
            int max_n = indexVec.size();
            // send node coefs from Phi[j] to bank
            // except for the root nodes, only wavelets are sent
            for (int i = 0; i < max_n; i++) {
                if (indexVec[i] < 0) continue; // nodes that are not in refOrb
                int csize = sizecoeffW;
                if (parindexVec[i] < 0) csize = sizecoeff;
                if (sizes > 0) { // fixed size
                    account.put_nodedata(j, indexVec[i], sizes, coeffVec[i]);
                } else {
                    account.put_nodedata(j, indexVec[i], csize, &(coeffVec[i][sizecoeff - csize]));
                }
            }
        }
        // Imaginary parts are considered as orbitals with an orbid shifted by N
        if (Phi[j]->hasImag()) {
            Phi[j]->imag().makeCoeffVector(coeffVec, indexVec, parindexVec, scalefac, max_ix, refTree);
            int max_n = indexVec.size();
            // send node coefs from Phi[j] to bank
            for (int i = 0; i < max_n; i++) {
                if (indexVec[i] < 0) continue; // nodes that are not in refOrb
                // NB: the identifier (indexVec[i]) must be shifted for not colliding with the nodes from the real part
                int csize = sizecoeffW;
                if (parindexVec[i] < 0) csize = sizecoeff;
                if (sizes > 0) { // fixed size
                    account.put_nodedata(j + N, indexVec[i], sizes, coeffVec[i]);
                } else {
                    account.put_nodedata(j + N, indexVec[i], csize, &(coeffVec[i][sizecoeff - csize]));
                }
            }
        }
    }
}

/** @brief Multiply all orbitals with a function
 *
 * @param Phi: orbitals to multiply
 * @param f  : function to multiply
 *
 * Computes the product of each orbital with a function
 * in parallel using a local representation.
 * Input trees are extended by one scale at most.
 */
CompFunctionVector multiply(CompFunctionVector &Phi, RepresentableFunction<3> &f, double prec, CompFunction<3> *Func, int nrefine, bool all) {

    int N = Phi.size();
    const int D = 3;
    bool serial = mpi::wrk_size == 1; // flag for serial/MPI switch

    // 1a) extend grid where f is large (around nuclei)
    // TODO: do it in save_nodes + refTree, only saving the extra nodes, without keeping them permanently. Or refine refTree?

    for (int i = 0; i < N; i++) {
        if (!mpi::my_func(i)) continue;
        int irefine = 0;
        while (Phi[i]->hasReal() and irefine < nrefine and refine_grid(Phi[i]->real(), f) > 0) irefine++;
        irefine = 0;
        while (Phi[i]->hasImag() and irefine < nrefine and refine_grid(Phi[i]->imag(), f) > 0) irefine++;
    }

    // 1b) make union tree without coefficients
    FunctionTree<D> refTree(*Phi.vecMRA);
    // refine_grid(refTree, f); //to test
    mpi::allreduce_Tree_noCoeff(refTree, Phi, mpi::comm_wrk);

    int kp1 = refTree.getKp1();
    int kp1_d = refTree.getKp1_d();
    int nCoefs = refTree.getTDim() * kp1_d;

    IntVector PsihasReIm = IntVector::Zero(2);
    for (int i = 0; i < N; i++) {
        if (!mpi::my_func(i)) continue;
        PsihasReIm[0] = (Phi[i]->hasReal()) ? 1 : 0;
        PsihasReIm[1] = (Phi[i]->hasImag()) ? 1 : 0;
    }
    mpi::allreduce_vector(PsihasReIm, mpi::comm_wrk);
    CompFunctionVector out(N);
    CompFunctionVector outtest(N);
    if (not PsihasReIm[0] and not PsihasReIm[1]) {
        return out; // do nothing
    }

    int Neff = N;
    if (PsihasReIm[1]) Neff = 2 * N; // Imag and Real treated independently. We always treat real part of Psi

    std::vector<double> scalefac_ref;
    std::vector<double *> coeffVec_ref; // not used!
    std::vector<int> indexVec_ref;      // serialIx of the nodes
    std::vector<int> parindexVec_ref;   // serialIx of the parent nodes
    std::vector<MWNode<D> *> refNodes;  // pointers to nodes
    int max_ix;
    // get a list of all nodes in union tree, identified by their serialIx indices
    refTree.makeCoeffVector(coeffVec_ref, indexVec_ref, parindexVec_ref, scalefac_ref, max_ix, refTree, &refNodes);
    int max_n = indexVec_ref.size();
    std::map<int, int> ix2n; // for a given serialIx, give index in vectors
    for (int nn = 0; nn < max_n; nn++) ix2n[indexVec_ref[nn]] = nn;

    // 2a) send own nodes to bank, identifying them through the serialIx of refTree
    BankAccount nodesPhi;        // to put the original nodes
    BankAccount nodesMultiplied; // to put the multiplied nodes

    // used for serial only:
    std::vector<std::vector<double *>> coeffVec(Neff);
    std::vector<std::vector<int>> indexVec(Neff);   // serialIx of the nodes
    std::map<int, std::vector<int>> node2orbVec;    // for each node index, gives a vector with the indices of the orbitals using this node
    std::vector<std::map<int, int>> orb2node(Neff); // for a given orbital and a given node, gives the node index in the
                                                    // orbital given the node index in the reference tree
    if (serial) {
        // make list of all coefficients (coeffVec), and their reference indices (indexVec)
        std::vector<int> parindexVec; // serialIx of the parent nodes
        std::vector<double> scalefac;
        for (int j = 0; j < N; j++) {
            // make vector with all coef pointers and their indices in the union grid
            if (Phi[j]->hasReal()) {
                Phi[j]->real().makeCoeffVector(coeffVec[j], indexVec[j], parindexVec, scalefac, max_ix, refTree);
                // make a map that gives j from indexVec
                int orb_node_ix = 0;
                for (int ix : indexVec[j]) {
                    orb2node[j][ix] = orb_node_ix++;
                    if (ix < 0) continue;
                    node2orbVec[ix].push_back(j);
                }
            }
            if (Phi[j]->hasImag()) {
                Phi[j]->imag().makeCoeffVector(coeffVec[j + N], indexVec[j + N], parindexVec, scalefac, max_ix, refTree);
                // make a map that gives j from indexVec
                int orb_node_ix = 0;
                for (int ix : indexVec[j + N]) {
                    orb2node[j + N][ix] = orb_node_ix++;
                    if (ix < 0) continue;
                    node2orbVec[ix].push_back(j + N);
                }
            }
        }
    } else {
        save_nodes(Phi, refTree, nodesPhi, nCoefs);
        mpi::barrier(mpi::comm_wrk); // required for now, as the blockdata functionality has no queue yet.
    }

    // 2b) save Func in bank and remove its coefficients
    if (Func != nullptr and !serial) {
        // put Func in local representation if not already done
        if (!Func->real().isLocal) { Func->real().saveNodesAndRmCoeff(); }
    }

    // 3) mutiply for each node
    std::vector<std::vector<double *>> coeffpVec(Neff); // to put pointers to the multiplied coefficient for each orbital in serial case
    std::vector<DoubleMatrix> multipliedCoeffVec;       // just to ensure that the data from multipliedCoeff is not deleted, since we point to it.
    std::vector<std::map<int, int>> ix2coef(Neff);      // to find the index in for example rotCoeffVec[] corresponding to a serialIx
    DoubleVector NODEP = DoubleVector::Zero(nCoefs);
    DoubleVector NODEF = DoubleVector::Zero(nCoefs);

    if (serial) {
#pragma omp parallel for schedule(dynamic)
        for (int n = 0; n < max_n; n++) {
            MWNode<D> node(*(refNodes[n]), false);
            int node_ix = indexVec_ref[n]; // SerialIx for this node in the reference tree

            // 3a) make values for f at this node
            // 3a1) get coordinates of quadrature points for this node
            Eigen::MatrixXd pts; // Eigen::Zero(D, nCoefs);
            double fval[nCoefs];
            Coord<D> r;
            double *originalCoef = nullptr;
            MWNode<3> *Fnode = nullptr;
            if (Func == nullptr) {
                node.getExpandedChildPts(pts); // TODO: use getPrimitiveChildPts (less cache).
                for (int j = 0; j < nCoefs; j++) {
                    for (int d = 0; d < D; d++) r[d] = pts(d, j); //*scaling_factor[d]?
                    fval[j] = f.evalf(r);
                }
            } else {
                Fnode = Func->real().findNode(node.getNodeIndex());
                if (Fnode == nullptr) {
                    node.getExpandedChildPts(pts); // TODO: use getPrimitiveChildPts (less cache).
                    for (int j = 0; j < nCoefs; j++) {
                        for (int d = 0; d < D; d++) r[d] = pts(d, j); //*scaling_factor[d]?
                        fval[j] = f.evalf(r);
                    }
                } else {
                    originalCoef = Fnode->getCoefs();
                    for (int j = 0; j < nCoefs; j++) fval[j] = originalCoef[j];
                    Fnode->attachCoefs(fval); // note that each thread has its own copy
                    Fnode->mwTransform(Reconstruction);
                    Fnode->cvTransform(Forward);
                }
            }
            DoubleMatrix multipliedCoeff(nCoefs, node2orbVec[node_ix].size());
            int i = 0;
            // 3b) fetch all orbitals at this node
            std::vector<int> orbjVec;            // to remember which orbital correspond to each orbVec.size();
            for (int j : node2orbVec[node_ix]) { // loop over indices of the orbitals using this node
                int orb_node_ix = orb2node[j][node_ix];
                orbjVec.push_back(j);
                for (int k = 0; k < nCoefs; k++) multipliedCoeff(k, i) = coeffVec[j][orb_node_ix][k];
                // 3c) transform to grid
                node.attachCoefs(&(multipliedCoeff(0, i)));
                node.mwTransform(Reconstruction);
                node.cvTransform(Forward);
                // 3d) multiply
                for (int k = 0; k < nCoefs; k++) multipliedCoeff(k, i) *= fval[k]; // replace by Matrix vector multiplication?
                // 3e) transform back to mw
                node.cvTransform(Backward);
                node.mwTransform(Compression);
                i++;
            }
            if (Func != nullptr and originalCoef != nullptr) {
                // restablish original values
                Fnode->attachCoefs(originalCoef);
            }

            // 3f) save multiplied nodes
            for (int i = 0; i < orbjVec.size(); i++) {
#pragma omp critical
                {
                    ix2coef[orbjVec[i]][node_ix] = coeffpVec[orbjVec[i]].size();
                    coeffpVec[orbjVec[i]].push_back(&(multipliedCoeff(0, i))); // list of coefficient pointers
                }
            }
#pragma omp critical
            {
                // this ensures that multipliedCoeff is not deleted, when getting out of scope
                multipliedCoeffVec.push_back(std::move(multipliedCoeff));
            }
            node.attachCoefs(nullptr); // to avoid deletion of valid multipliedCoeff by destructor
        }
    } else {
        // MPI
        int count1 = 0;
        int count2 = 0;
        TaskManager tasks(max_n);
        for (int nn = 0; nn < max_n; nn++) {
            int n = tasks.next_task();
            if (n < 0) break;
            MWNode<D> node(*(refNodes[n]), false);
            // 3a) make values for f
            // 3a1) get coordinates of quadrature points for this node
            Eigen::MatrixXd pts;           // Eigen::Zero(D, nCoefs);
            node.getExpandedChildPts(pts); // TODO: use getPrimitiveChildPts (less cache).
            double fval[nCoefs];
            Coord<D> r;
            MWNode<D> Fnode(*(refNodes[n]), false);
            if (Func == nullptr) {
                for (int j = 0; j < nCoefs; j++) {
                    for (int d = 0; d < D; d++) r[d] = pts(d, j); //*scaling_factor[d]?
                    fval[j] = f.evalf(r);
                }
            } else {
                int nIdx = Func->real().getIx(node.getNodeIndex());
                count1++;
                if (nIdx < 0) {
                    // use the function f instead of Func
                    count2++;
                    for (int j = 0; j < nCoefs; j++) {
                        for (int d = 0; d < D; d++) r[d] = pts(d, j);
                        fval[j] = f.evalf(r);
                    }
                } else {
                    Func->real().getNodeCoeff(nIdx, fval); // fetch coef from Bank
                    Fnode.attachCoefs(fval);
                    Fnode.mwTransform(Reconstruction);
                    Fnode.cvTransform(Forward);
                }
            }

            // 3b) fetch all orbitals at this node
            DoubleMatrix coeffBlock(nCoefs, Neff); // largest possible used size
            std::vector<int> orbjVec;
            nodesPhi.get_nodeblock(indexVec_ref[n], coeffBlock.data(), orbjVec);
            coeffBlock.conservativeResize(Eigen::NoChange, orbjVec.size()); // keep only used part
            DoubleMatrix MultipliedCoeff(nCoefs, orbjVec.size());
            // 3c) transform to grid
            for (int j = 0; j < orbjVec.size(); j++) { // TODO: transform all j at once ?
                // TODO: select only nodes that are end nodes?
                node.attachCoefs(coeffBlock.col(j).data());
                node.mwTransform(Reconstruction);
                node.cvTransform(Forward);
                // 3d) multiply
                double *coefs = node.getCoefs();
                for (int i = 0; i < nCoefs; i++) coefs[i] *= fval[i];
                // 3e) transform back to mw
                node.cvTransform(Backward);
                node.mwTransform(Compression);
                // 3f) save multiplied nodes
                nodesMultiplied.put_nodedata(orbjVec[j], indexVec_ref[n] + max_ix, nCoefs, coefs);
            }
            node.attachCoefs(nullptr);  // to avoid deletion of valid multipliedCoeff by destructor
            Fnode.attachCoefs(nullptr); // to avoid deletion of valid multipliedCoeff by destructor
        }
        mrcpp::mpi::barrier(mrcpp::mpi::comm_wrk); // wait until everything is stored before fetching!
    }

    // 5) reconstruct trees using multiplied nodes.

    // only serial case can use OMP, because MPI cannot be used by threads
    if (serial) {
        // OMP parallelized, but does not scale well, because the total memory bandwidth is a bottleneck. (the main
        // operation is writing the coefficient into the tree)

#pragma omp parallel for schedule(static)
        for (int j = 0; j < Neff; j++) {
            if (j < N) {
                if (Phi[j]->hasReal()) {
                    out[j]->alloc(0);
                    out[j]->real().clear();
                    out[j]->real().makeTreefromCoeff(refTree, coeffpVec[j], ix2coef[j], -1.0, "copy");
                    // 6) reconstruct trees from end nodes
                    out[j]->real().mwTransform(BottomUp);
                    out[j]->real().calcSquareNorm();
                }
            } else {
                if (Phi[j % N]->hasImag()) {
                    out[j % N]->alloc(0);
                    out[j % N]->imag().clear();
                    out[j % N]->imag().makeTreefromCoeff(refTree, coeffpVec[j], ix2coef[j], -1.0, "copy");
                    out[j]->imag().mwTransform(BottomUp);
                    out[j]->imag().calcSquareNorm();
                }
            }
        }
    } else {
        for (int j = 0; j < Neff; j++) {
            if (not mpi::my_func(j % N) and not all) continue;
            // traverse possible nodes, and stop descending when norm is zero (leaf in out[j])
            std::vector<double *> coeffpVec; //
            std::map<int, int> ix2coef;      // to find the index in coeffVec[] corresponding to a serialIx in refTree
            int ix = 0;
            std::vector<double *> pointerstodelete; // list of temporary arrays to clean up

            for (int ibank = 0; ibank < mpi::bank_size; ibank++) {
                std::vector<int> nodeidVec;
                double *dataVec; // will be allocated by bank
                nodesMultiplied.get_orbblock(j, dataVec, nodeidVec, ibank);
                if (nodeidVec.size() > 0) pointerstodelete.push_back(dataVec);
                int shift = 0;
                for (int n = 0; n < nodeidVec.size(); n++) {
                    assert(nodeidVec[n] - max_ix >= 0);                // unmultiplied nodes have been deleted
                    assert(ix2coef.count(nodeidVec[n] - max_ix) == 0); // each nodeid treated once
                    ix2coef[nodeidVec[n] - max_ix] = ix++;
                    coeffpVec.push_back(&dataVec[shift]); // list of coeff pointers
                    shift += nCoefs;
                }
            }
            if (j < N) {
                if (Phi[j]->hasReal()) {
                    out[j]->alloc(0);
                    out[j]->real().clear();
                    out[j]->real().makeTreefromCoeff(refTree, coeffpVec, ix2coef, -1.0, "copy");
                    // 6) reconstruct trees from end nodes
                    out[j]->real().mwTransform(BottomUp);
                    out[j]->real().calcSquareNorm();
                    out[j]->real().resetEndNodeTable();
                    // out[j]->real().crop(prec, 1.0, false); //bad convergence if out is cropped
                    if (nrefine > 0) Phi[j]->real().crop(prec, 1.0, false); // restablishes original Phi
                }
            } else {
                if (Phi[j % N]->hasImag()) {
                    out[j % N]->alloc(0);
                    out[j % N]->imag().clear();
                    out[j % N]->imag().makeTreefromCoeff(refTree, coeffpVec, ix2coef, -1.0, "copy");
                    out[j % N]->imag().mwTransform(BottomUp);
                    out[j % N]->imag().calcSquareNorm();
                    // out[j % N]->imag().crop(prec, 1.0, false);
                    if (nrefine > 0) Phi[j % N]->imag().crop(prec, 1.0, false);
                }
            }

            for (double *p : pointerstodelete) delete[] p;
            pointerstodelete.clear();
        }
    }
    return out;
}

void SetdefaultMRA(MultiResolutionAnalysis<3> *MRA) {
    defaultCompMRA<3> = MRA;
}

ComplexVector dot(CompFunctionVector &Bra, CompFunctionVector &Ket) {
    int N = Bra.size();
    ComplexVector result = ComplexVector::Zero(N);
    for (int i = 0; i < N; i++) {
        // The bra is sent to the owner of the ket
        if (my_func(Bra[i]) != my_func(Ket[i])) { MSG_ABORT("same indices should have same ownership"); }
        result[i] = dot(*Bra[i], *Ket[i]);
        if (not mrcpp::mpi::my_func(i)) Bra[i]->free();
    }
    mrcpp::mpi::allreduce_vector(result, mrcpp::mpi::comm_wrk);
    return result;
}

/** @brief Compute Löwdin orthonormalization matrix
 *
 * @param Phi: orbitals to orthonomalize
 *
 * Computes the inverse square root of the orbital overlap matrix S^(-1/2)
 */
ComplexMatrix calc_lowdin_matrix(CompFunctionVector &Phi) {
    ComplexMatrix S_tilde = calc_overlap_matrix(Phi);
    ComplexMatrix S_m12 = math_utils::hermitian_matrix_pow(S_tilde, -1.0 / 2.0);
    return S_m12;
}

/** @brief Orbital transformation out_j = sum_i inp_i*U_ij
 *
 * NOTE: OrbitalVector is considered a ROW vector, so rotation
 *       means matrix multiplication from the right
 *
 * MPI: Rank distribution of output vector is the same as input vector
 *
 */
ComplexMatrix calc_overlap_matrix(CompFunctionVector &BraKet) {
    // NB: must be spinseparated at this point!

    int N = BraKet.size();
    ComplexMatrix S = ComplexMatrix::Zero(N, N);
    DoubleMatrix Sreal = DoubleMatrix::Zero(2 * N, 2 * N); // same as S, but stored as 4 blocks, rr,ri,ir,ii
    MultiResolutionAnalysis<3> *mra = BraKet.vecMRA;

    // 1) make union tree without coefficients
    mrcpp::FunctionTree<3> refTree(*mra);
    mpi::allreduce_Tree_noCoeff(refTree, BraKet, mpi::comm_wrk);

    int sizecoeff = (1 << refTree.getDim()) * refTree.getKp1_d();
    int sizecoeffW = ((1 << refTree.getDim()) - 1) * refTree.getKp1_d();

    // get a list of all nodes in union grid, as defined by their indices
    std::vector<double> scalefac;
    std::vector<double *> coeffVec_ref;
    std::vector<int> indexVec_ref;    // serialIx of the nodes
    std::vector<int> parindexVec_ref; // serialIx of the parent nodes
    int max_ix;                       // largest index value (not used here)

    refTree.makeCoeffVector(coeffVec_ref, indexVec_ref, parindexVec_ref, scalefac, max_ix, refTree);
    int max_n = indexVec_ref.size();

    // only used for serial case:
    std::vector<std::vector<double *>> coeffVec(2 * N);
    std::map<int, std::vector<int>> node2orbVec;     // for each node index, gives a vector with the indices of the orbitals using this node
    std::vector<std::map<int, int>> orb2node(2 * N); // for a given orbital and a given node, gives the node index in
                                                     // the orbital given the node index in the reference tree

    bool serial = mrcpp::mpi::wrk_size == 1; // flag for serial/MPI switch
    mrcpp::BankAccount nodesBraKet;

    // In the serial case we store the coeff pointers in coeffVec. In the mpi case the coeff are stored in the bank
    if (serial) {
        // 2) make list of all coefficients, and their reference indices
        // for different orbitals, indexVec will give the same index for the same node in space
        std::vector<int> parindexVec; // serialIx of the parent nodes
        std::vector<int> indexVec;    // serialIx of the nodes
        for (int j = 0; j < N; j++) {
            // make vector with all coef pointers and their indices in the union grid
            if (BraKet[j]->hasReal()) {
                BraKet[j]->real().makeCoeffVector(coeffVec[j], indexVec, parindexVec, scalefac, max_ix, refTree);
                // make a map that gives j from indexVec
                int orb_node_ix = 0;
                for (int ix : indexVec) {
                    orb2node[j][ix] = orb_node_ix++;
                    if (ix < 0) continue;
                    node2orbVec[ix].push_back(j);
                }
            }
            if (BraKet[j]->hasImag()) {
                BraKet[j]->imag().makeCoeffVector(coeffVec[j + N], indexVec, parindexVec, scalefac, max_ix, refTree);
                // make a map that gives j from indexVec
                int orb_node_ix = 0;
                for (int ix : indexVec) {
                    orb2node[j + N][ix] = orb_node_ix++;
                    if (ix < 0) continue;
                    node2orbVec[ix].push_back(j + N);
                }
            }
        }
    } else { // MPI case
        // 2) send own nodes to bank, identifying them through the serialIx of refTree
        save_nodes(BraKet, refTree, nodesBraKet);
        mrcpp::mpi::barrier(mrcpp::mpi::comm_wrk); // wait until everything is stored before fetching!
    }

    // 3) make dot product for all the nodes and accumulate into S

    int ibank = 0;
#pragma omp parallel for schedule(dynamic) if (serial)
    for (int n = 0; n < max_n; n++) {
        if (n % mrcpp::mpi::wrk_size != mrcpp::mpi::wrk_rank) continue;
        int csize;
        int node_ix = indexVec_ref[n]; // SerialIx for this node in the reference tree
        std::vector<int> orbVec;       // identifies which orbitals use this node
        if (serial and node2orbVec[node_ix].size() <= 0) continue;
        if (parindexVec_ref[n] < 0)
            csize = sizecoeff;
        else
            csize = sizecoeffW;

        // In the serial case we copy the coeff coeffBlock. In the mpi case coeffBlock is provided by the bank
        if (serial) {
            int shift = sizecoeff - sizecoeffW; // to copy only wavelet part
            if (parindexVec_ref[n] < 0) shift = 0;
            DoubleMatrix coeffBlock(csize, node2orbVec[node_ix].size());
            for (int j : node2orbVec[node_ix]) { // loop over indices of the orbitals using this node
                int orb_node_ix = orb2node[j][node_ix];
                for (int k = 0; k < csize; k++) coeffBlock(k, orbVec.size()) = coeffVec[j][orb_node_ix][k + shift];
                orbVec.push_back(j);
            }
            if (orbVec.size() > 0) {
                DoubleMatrix S_temp(orbVec.size(), orbVec.size());
                S_temp.noalias() = coeffBlock.transpose() * coeffBlock;
                for (int i = 0; i < orbVec.size(); i++) {
                    for (int j = 0; j < orbVec.size(); j++) {
                        if (BraKet[orbVec[i] % N]->data.n1[0] != BraKet[orbVec[j] % N]->data.n1[0] and
                            BraKet[orbVec[i] % N]->data.n1[0] != 0 and
                            BraKet[orbVec[j] % N]->data.n1[0] != 0)
                            continue;
                        double &Srealij = Sreal(orbVec[i], orbVec[j]);
                        double &Stempij = S_temp(i, j);
#pragma omp atomic
                        Srealij += Stempij;
                    }
                }
            }
        } else { // MPI case
            DoubleMatrix coeffBlock(csize, 2 * N);
            nodesBraKet.get_nodeblock(indexVec_ref[n], coeffBlock.data(), orbVec);

            if (orbVec.size() > 0) {
                DoubleMatrix S_temp(orbVec.size(), orbVec.size());
                coeffBlock.conservativeResize(Eigen::NoChange, orbVec.size());
                S_temp.noalias() = coeffBlock.transpose() * coeffBlock;
                for (int i = 0; i < orbVec.size(); i++) {
                    for (int j = 0; j < orbVec.size(); j++) {
                        if (BraKet[orbVec[i] % N]->data.n1[0] != BraKet[orbVec[j] % N]->data.n1[0] and
                            BraKet[orbVec[i] % N]->data.n1[0] != 0 and
                            BraKet[orbVec[j] % N]->data.n1[0] != 0)
                            continue;
                        Sreal(orbVec[i], orbVec[j]) += S_temp(i, j);
                    }
                }
            }
        }
    }
    IntVector conjMat = IntVector::Zero(N);
    for (int i = 0; i < N; i++) {
        if (!mrcpp::mpi::my_func(BraKet[i])) continue;
        conjMat[i] = (BraKet[i]->conjugate()) ? -1 : 1;
    }
    mrcpp::mpi::allreduce_vector(conjMat, mrcpp::mpi::comm_wrk);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j <= i; j++) {
            S.real()(i, j) = Sreal(i, j) + conjMat[i] * conjMat[j] * Sreal(i + N, j + N);
            S.imag()(i, j) = conjMat[j] * Sreal(i, j + N) - conjMat[i] * Sreal(i + N, j);
            if (i != j) S(j, i) = std::conj(S(i, j)); // ensure exact symmetri
        }
    }

    // Assumes linearity: result is sum of all nodes contributions
    mrcpp::mpi::allreduce_matrix(S, mrcpp::mpi::comm_wrk);

    return S;
}

/** @brief Compute the overlap matrix S_ij = <bra_i|ket_j>
 *
 */
ComplexMatrix calc_overlap_matrix(CompFunctionVector &Bra, CompFunctionVector &Ket) {
    mrcpp::mpi::barrier(mrcpp::mpi::comm_wrk); // for consistent timings

    MultiResolutionAnalysis<3> *mra = Bra.vecMRA;

    int N = Bra.size();
    int M = Ket.size();
    ComplexMatrix S = ComplexMatrix::Zero(N, M);
    DoubleMatrix Sreal = DoubleMatrix::Zero(2 * N, 2 * M); // same as S, but stored as 4 blocks, rr,ri,ir,ii

    // 1) make union tree without coefficients for Bra (supposed smallest)
    mrcpp::FunctionTree<3> refTree(*mra);
    mrcpp::mpi::allreduce_Tree_noCoeff(refTree, Bra, mpi::comm_wrk);
    // note that Ket is not part of union grid: if a node is in ket but not in Bra, the dot product is zero.

    int sizecoeff = (1 << refTree.getDim()) * refTree.getKp1_d();
    int sizecoeffW = ((1 << refTree.getDim()) - 1) * refTree.getKp1_d();

    // get a list of all nodes in union grid, as defined by their indices
    std::vector<double *> coeffVec_ref;
    std::vector<int> indexVec_ref;    // serialIx of the nodes
    std::vector<int> parindexVec_ref; // serialIx of the parent nodes
    std::vector<double> scalefac;
    int max_ix;

    refTree.makeCoeffVector(coeffVec_ref, indexVec_ref, parindexVec_ref, scalefac, max_ix, refTree);
    int max_n = indexVec_ref.size();
    max_ix++;

    bool serial = mrcpp::mpi::wrk_size == 1; // flag for serial/MPI switch

    // only used for serial case:
    std::vector<std::vector<double *>> coeffVecBra(2 * N);
    std::map<int, std::vector<int>> node2orbVecBra;     // for each node index, gives a vector with the indices of the orbitals using this node
    std::vector<std::map<int, int>> orb2nodeBra(2 * N); // for a given orbital and a given node, gives the node index in
                                                        // the orbital given the node index in the reference tree
    std::vector<std::vector<double *>> coeffVecKet(2 * M);
    std::map<int, std::vector<int>> node2orbVecKet;     // for each node index, gives a vector with the indices of the orbitals using this node
    std::vector<std::map<int, int>> orb2nodeKet(2 * M); // for a given orbital and a given node, gives the node index in
                                                        // the orbital given the node index in the reference tree
    mrcpp::BankAccount nodesBra;
    mrcpp::BankAccount nodesKet;

    // In the serial case we store the coeff pointers in coeffVec. In the mpi case the coeff are stored in the bank
    if (serial) {
        // 2) make list of all coefficients, and their reference indices
        // for different orbitals, indexVec will give the same index for the same node in space
        // TODO? : do not copy coefficients, but use directly the pointers
        // could OMP parallelize, but is fast anyway
        std::vector<int> parindexVec; // serialIx of the parent nodes
        std::vector<int> indexVec;    // serialIx of the nodes
        for (int j = 0; j < N; j++) {
            // make vector with all coef pointers and their indices in the union grid
            if (Bra[j]->hasReal()) {
                Bra[j]->real().makeCoeffVector(coeffVecBra[j], indexVec, parindexVec, scalefac, max_ix, refTree);
                // make a map that gives j from indexVec
                int orb_node_ix = 0;
                for (int ix : indexVec) {
                    orb2nodeBra[j][ix] = orb_node_ix++;
                    if (ix < 0) continue;
                    node2orbVecBra[ix].push_back(j);
                }
            }
            if (Bra[j]->hasImag()) {
                Bra[j]->imag().makeCoeffVector(coeffVecBra[j + N], indexVec, parindexVec, scalefac, max_ix, refTree);
                // make a map that gives j from indexVec
                int orb_node_ix = 0;
                for (int ix : indexVec) {
                    orb2nodeBra[j + N][ix] = orb_node_ix++;
                    if (ix < 0) continue;
                    node2orbVecBra[ix].push_back(j + N);
                }
            }
        }
        for (int j = 0; j < M; j++) {
            if (Ket[j]->hasReal()) {
                Ket[j]->real().makeCoeffVector(coeffVecKet[j], indexVec, parindexVec, scalefac, max_ix, refTree);
                // make a map that gives j from indexVec
                int orb_node_ix = 0;
                for (int ix : indexVec) {
                    orb2nodeKet[j][ix] = orb_node_ix++;
                    if (ix < 0) continue;
                    node2orbVecKet[ix].push_back(j);
                }
            }
            if (Ket[j]->hasImag()) {
                Ket[j]->imag().makeCoeffVector(coeffVecKet[j + M], indexVec, parindexVec, scalefac, max_ix, refTree);
                // make a map that gives j from indexVec
                int orb_node_ix = 0;
                for (int ix : indexVec) {
                    orb2nodeKet[j + M][ix] = orb_node_ix++;
                    if (ix < 0) continue;
                    node2orbVecKet[ix].push_back(j + M);
                }
            }
        }

    } else { // MPI case
        // 2) send own nodes to bank, identifying them through the serialIx of refTree
        save_nodes(Bra, refTree, nodesBra);
        save_nodes(Ket, refTree, nodesKet);
        mrcpp::mpi::barrier(mrcpp::mpi::comm_wrk); // wait until everything is stored before fetching!
    }

    // 3) make dot product for all the nodes and accumulate into S
    int totsiz = 0;
    int totget = 0;
    int mxtotsiz = 0;
    int ibank = 0;
    //For some unknown reason the h2_mag_lda test sometimes fails when schedule(dynamic) is chosen
#pragma omp parallel for schedule(static) if (serial)
    for (int n = 0; n < max_n; n++) {
        if (n % mrcpp::mpi::wrk_size != mrcpp::mpi::wrk_rank) continue;
        int csize;
        std::vector<int> orbVecBra; // identifies which Bra orbitals use this node
        std::vector<int> orbVecKet; // identifies which Ket orbitals use this node
        if (parindexVec_ref[n] < 0)
            csize = sizecoeff;
        else
            csize = sizecoeffW;
        if (serial) {
            int node_ix = indexVec_ref[n];      // SerialIx for this node in the reference tree
            int shift = sizecoeff - sizecoeffW; // to copy only wavelet part
            DoubleMatrix coeffBlockBra(csize, node2orbVecBra[node_ix].size());
            DoubleMatrix coeffBlockKet(csize, node2orbVecKet[node_ix].size());
            if (parindexVec_ref[n] < 0) shift = 0;

            for (int j : node2orbVecBra[node_ix]) { // loop over indices of the orbitals using this node
                int orb_node_ix = orb2nodeBra[j][node_ix];
                for (int k = 0; k < csize; k++) coeffBlockBra(k, orbVecBra.size()) = coeffVecBra[j][orb_node_ix][k + shift];
                orbVecBra.push_back(j);
            }
            for (int j : node2orbVecKet[node_ix]) { // loop over indices of the orbitals using this node
                int orb_node_ix = orb2nodeKet[j][node_ix];
                for (int k = 0; k < csize; k++) coeffBlockKet(k, orbVecKet.size()) = coeffVecKet[j][orb_node_ix][k + shift];
                orbVecKet.push_back(j);
            }

            if (orbVecBra.size() > 0 and orbVecKet.size() > 0) {
                DoubleMatrix S_temp(orbVecBra.size(), orbVecKet.size());
                S_temp.noalias() = coeffBlockBra.transpose() * coeffBlockKet;
                for (int i = 0; i < orbVecBra.size(); i++) {
                    for (int j = 0; j < orbVecKet.size(); j++) {
                        if (Bra[orbVecBra[i] % N]->data.n1[0] != Ket[orbVecKet[j] % M]->data.n1[0] and
                            Bra[orbVecBra[i] % N]->data.n1[0] != 0 and
                            Ket[orbVecKet[j] % M]->data.n1[0] != 0)
                            continue;
                        // must ensure that threads are not competing
                        double &Srealij = Sreal(orbVecBra[i], orbVecKet[j]);
                        double &Stempij = S_temp(i, j);
#pragma omp atomic
                        Srealij += Stempij;
                    }
                }
            }
        } else {

            DoubleMatrix coeffBlockBra(csize, 2 * N);
            DoubleMatrix coeffBlockKet(csize, 2 * M);
            nodesBra.get_nodeblock(indexVec_ref[n], coeffBlockBra.data(), orbVecBra); // get Bra parts
            nodesKet.get_nodeblock(indexVec_ref[n], coeffBlockKet.data(), orbVecKet); // get Ket parts
            totsiz += orbVecBra.size() * orbVecKet.size();
            mxtotsiz += N * M;
            totget += orbVecBra.size() + orbVecKet.size();
            if (orbVecBra.size() > 0 and orbVecKet.size() > 0) {
                DoubleMatrix S_temp(orbVecBra.size(), orbVecKet.size());
                coeffBlockBra.conservativeResize(Eigen::NoChange, orbVecBra.size());
                coeffBlockKet.conservativeResize(Eigen::NoChange, orbVecKet.size());
                S_temp.noalias() = coeffBlockBra.transpose() * coeffBlockKet;
                for (int i = 0; i < orbVecBra.size(); i++) {
                    for (int j = 0; j < orbVecKet.size(); j++) {
                        if (Bra[orbVecBra[i] % N]->data.n1[0] != Ket[orbVecKet[j] % M]->data.n1[0] and
                            Bra[orbVecBra[i] % N]->data.n1[0] != 0 and
                            Ket[orbVecKet[j] % M]->data.n1[0] != 0)
                            continue;
                        Sreal(orbVecBra[i], orbVecKet[j]) += S_temp(i, j);
                    }
                }
            }
        }
    }

    IntVector conjMatBra = IntVector::Zero(N);
    for (int i = 0; i < N; i++) {
        if (!mrcpp::mpi::my_func(Bra[i])) continue;
        conjMatBra[i] = (Bra[i]->conjugate()) ? -1 : 1;
    }
    mrcpp::mpi::allreduce_vector(conjMatBra, mrcpp::mpi::comm_wrk);
    IntVector conjMatKet = IntVector::Zero(M);
    for (int i = 0; i < M; i++) {
        if (!mrcpp::mpi::my_func(Ket[i])) continue;
        conjMatKet[i] = (Ket[i]->conjugate()) ? -1 : 1;
    }
    mrcpp::mpi::allreduce_vector(conjMatKet, mrcpp::mpi::comm_wrk);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            S.real()(i, j) = Sreal(i, j) + conjMatBra[i] * conjMatKet[j] * Sreal(i + N, j + M);
            S.imag()(i, j) = conjMatKet[j] * Sreal(i, j + M) - conjMatBra[i] * Sreal(i + N, j);
        }
    }

    // 4) collect results from all MPI. Linearity: result is sum of all node contributions

    mrcpp::mpi::allreduce_matrix(S, mrcpp::mpi::comm_wrk);

    return S;
}

/** @brief Compute the overlap matrix of the absolute value of the functions S_ij = <|bra_i|||ket_j|>
 *
 */
DoubleMatrix calc_norm_overlap_matrix(CompFunctionVector &BraKet) {
    int N = BraKet.size();
    DoubleMatrix S = DoubleMatrix::Zero(N, N);
    DoubleMatrix Sreal = DoubleMatrix::Zero(2 * N, 2 * N); // same as S, but stored as 4 blocks, rr,ri,ir,ii
    MultiResolutionAnalysis<3> *mra = BraKet.vecMRA;

    // 1) make union tree without coefficients
    mrcpp::FunctionTree<3> refTree(*mra);
    mrcpp::mpi::allreduce_Tree_noCoeff(refTree, BraKet, mpi::comm_wrk);

    int sizecoeff = (1 << refTree.getDim()) * refTree.getKp1_d();
    int sizecoeffW = ((1 << refTree.getDim()) - 1) * refTree.getKp1_d();

    // get a list of all nodes in union grid, as defined by their indices
    std::vector<double> scalefac;
    std::vector<double *> coeffVec_ref;
    std::vector<int> indexVec_ref;    // serialIx of the nodes
    std::vector<int> parindexVec_ref; // serialIx of the parent nodes
    int max_ix;                       // largest index value (not used here)

    refTree.makeCoeffVector(coeffVec_ref, indexVec_ref, parindexVec_ref, scalefac, max_ix, refTree);
    int max_n = indexVec_ref.size();

    // only used for serial case:
    std::vector<std::vector<double *>> coeffVec(2 * N);
    std::map<int, std::vector<int>> node2orbVec;     // for each node index, gives a vector with the indices of the orbitals using this node
    std::vector<std::map<int, int>> orb2node(2 * N); // for a given orbital and a given node, gives the node index in
                                                     // the orbital given the node index in the reference tree

    bool serial = mrcpp::mpi::wrk_size == 1; // flag for serial/MPI switch
    mrcpp::BankAccount nodesBraKet;

    // In the serial case we store the coeff pointers in coeffVec. In the mpi case the coeff are stored in the bank
    if (serial) {
        // 2) make list of all coefficients, and their reference indices
        // for different orbitals, indexVec will give the same index for the same node in space
        std::vector<int> parindexVec; // serialIx of the parent nodes
        std::vector<int> indexVec;    // serialIx of the nodes
        for (int j = 0; j < N; j++) {
            // make vector with all coef pointers and their indices in the union grid
            if (BraKet[j]->hasReal()) {
                BraKet[j]->real().makeCoeffVector(coeffVec[j], indexVec, parindexVec, scalefac, max_ix, refTree);
                // make a map that gives j from indexVec
                int orb_node_ix = 0;
                for (int ix : indexVec) {
                    orb2node[j][ix] = orb_node_ix++;
                    if (ix < 0) continue;
                    node2orbVec[ix].push_back(j);
                }
            }
            if (BraKet[j]->hasImag()) {
                BraKet[j]->imag().makeCoeffVector(coeffVec[j + N], indexVec, parindexVec, scalefac, max_ix, refTree);
                // make a map that gives j from indexVec
                int orb_node_ix = 0;
                for (int ix : indexVec) {
                    orb2node[j + N][ix] = orb_node_ix++;
                    if (ix < 0) continue;
                    node2orbVec[ix].push_back(j + N);
                }
            }
        }
    } else { // MPI case
        // 2) send own nodes to bank, identifying them through the serialIx of refTree
        save_nodes(BraKet, refTree, nodesBraKet);
        mrcpp::mpi::barrier(mrcpp::mpi::comm_wrk); // wait until everything is stored before fetching!
    }

    // 3) make dot product for all the nodes and accumulate into S

    int ibank = 0;
#pragma omp parallel for schedule(dynamic) if (serial)
    for (int n = 0; n < max_n; n++) {
        if (n % mrcpp::mpi::wrk_size != mrcpp::mpi::wrk_rank) continue;
        int csize;
        int node_ix = indexVec_ref[n]; // SerialIx for this node in the reference tree
        std::vector<int> orbVec;       // identifies which orbitals use this node
        if (serial and node2orbVec[node_ix].size() <= 0) continue;
        if (parindexVec_ref[n] < 0)
            csize = sizecoeff;
        else
            csize = sizecoeffW;
        // In the serial case we copy the coeff coeffBlock. In the mpi case coeffBlock is provided by the bank
        if (serial) {
            int shift = sizecoeff - sizecoeffW; // to copy only wavelet part
            if (parindexVec_ref[n] < 0) shift = 0;
            DoubleMatrix coeffBlock(csize, node2orbVec[node_ix].size());
            for (int j : node2orbVec[node_ix]) { // loop over indices of the orbitals using this node
                int orb_node_ix = orb2node[j][node_ix];
                for (int k = 0; k < csize; k++) coeffBlock(k, orbVec.size()) = coeffVec[j][orb_node_ix][k + shift];
                orbVec.push_back(j);
            }
            if (orbVec.size() > 0) {
                DoubleMatrix S_temp(orbVec.size(), orbVec.size());
                coeffBlock = coeffBlock.cwiseAbs();
                S_temp.noalias() = coeffBlock.transpose() * coeffBlock;
                for (int i = 0; i < orbVec.size(); i++) {
                    for (int j = 0; j < orbVec.size(); j++) {
                        if (BraKet[orbVec[i] % N]->data.n1[0] != BraKet[orbVec[j] % N]->data.n1[0] and
                            BraKet[orbVec[i] % N]->data.n1[0] != 0 and
                            BraKet[orbVec[j] % N]->data.n1[0]!= 0)
                            continue;
                        double &Srealij = Sreal(orbVec[i], orbVec[j]);
                        double &Stempij = S_temp(i, j);
#pragma omp atomic
                        Srealij += Stempij;
                    }
                }
            }
        } else { // MPI case
            DoubleMatrix coeffBlock(csize, 2 * N);
            nodesBraKet.get_nodeblock(indexVec_ref[n], coeffBlock.data(), orbVec);

            if (orbVec.size() > 0) {
                DoubleMatrix S_temp(orbVec.size(), orbVec.size());
                coeffBlock.conservativeResize(Eigen::NoChange, orbVec.size());
                coeffBlock = coeffBlock.cwiseAbs();
                S_temp.noalias() = coeffBlock.transpose() * coeffBlock;
                for (int i = 0; i < orbVec.size(); i++) {
                    for (int j = 0; j < orbVec.size(); j++) {
                        if (BraKet[orbVec[i] % N]->data.n1[0] != BraKet[orbVec[j] % N]->data.n1[0] and
                            BraKet[orbVec[i] % N]->data.n1[0] != 0 and
                            BraKet[orbVec[j] % N]->data.n1[0]!= 0)
                            continue;
                        Sreal(orbVec[i], orbVec[j]) += S_temp(i, j);
                    }
                }
            }
        }
    }

    IntVector conjMat = IntVector::Zero(N);
    for (int i = 0; i < N; i++) {
        if (!mrcpp::mpi::my_func(i)) continue;
        conjMat[i] = (BraKet[i]->conjugate()) ? -1 : 1;
    }
    mrcpp::mpi::allreduce_vector(conjMat, mrcpp::mpi::comm_wrk);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j <= i; j++) {
            S(i, j) = Sreal(i, j) + conjMat[i] * conjMat[j] * Sreal(i + N, j + N) + conjMat[j] * Sreal(i, j + N) - conjMat[i] * Sreal(i + N, j);
            S(j, i) = S(i, j);
        }
    }

    // Assumes linearity: result is sum of all nodes contributions
    mrcpp::mpi::allreduce_matrix(S, mrcpp::mpi::comm_wrk);
    return S;
}

/** @brief Orthogonalize the functions in Bra against all orbitals in Ket
 *
 */
void orthogonalize(double prec, CompFunctionVector &Bra, CompFunctionVector &Ket) {
    // TODO: generalize for cases where Ket functions are not orthogonal to each other?
    ComplexMatrix S = calc_overlap_matrix(Bra, Ket);
    int N = Bra.size();
    int M = Ket.size();
    DoubleVector Ketnorms = DoubleVector::Zero(M);
    for (int i = 0; i < M; i++) {
        if (mpi::my_func(Ket[i])) Ketnorms(i)  = Ket[i]->squaredNorm();
    }
    mrcpp::mpi::allreduce_vector(Ketnorms, mrcpp::mpi::comm_wrk);
    ComplexMatrix rmat =  ComplexMatrix::Zero(M, N);
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < M; i++) {
            rmat(i,j) = 0.0 - S.conjugate()(j,i)/Ketnorms(i);
        }
    }
    CompFunctionVector rotatedKet(N);
    rotate(Ket, rmat, rotatedKet, prec / M);
    for (int j = 0; j < N; j++) {
        if(my_func(Bra[j]))Bra[j]->add(1.0,*rotatedKet[j]);
    }
}
template ComplexDouble dot(CompFunction<3> bra, CompFunction<3> ket) ;

} // namespace mrcpp
