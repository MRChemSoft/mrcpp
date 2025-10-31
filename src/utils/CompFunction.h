/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2021 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
 *
 * This file is part of MRCPP.
 *
 * MRCPP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRCPP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRCPP.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRCPP, see:
 * <https://mrcpp.readthedocs.io/>
 */

#pragma once
/**
 * @file
 * @brief Composite multicomponent function types (real/complex) on MRCPP multiresolution trees.
 *
 * This header defines:
 * - @ref CompFunctionData : POD metadata describing a multicomponent function.
 * - @ref TreePtr          : Small owning handle to up to four component trees (real/complex),
 *                           with optional MPI shared-memory backing.
 * - @ref CompFunction     : A high-level wrapper that owns/addresses component trees,
 *                           provides algebra (add/multiply/dot), projection, scaling,
 *                           norms, and utilities.
 * - Helpers for deep copies, linear combinations, products, projections, and orthogonalization.
 * - @ref CompFunctionVector : Convenience container for 3D functions with utilities for
 *                             rotations and overlap matrices.
 *
 * Components are stored as MRCPP @ref FunctionTree "FunctionTree<D, T>" instances.
 * Both real (`double`) and complex (`ComplexDouble`) representations are supported.
 *
 * Parallel notes:
 * - If built with MPI and `is_shared == true`, @ref TreePtr can allocate backing storage
 *   in an MPI shared-memory window (per @ref mpi::comm_share) to reduce duplication.
 * - Distribution utilities (e.g., @ref CompFunctionVector::distribute) use the runtime in
 *   `mpi_utils.h`.
 */

#include "mpi_utils.h"
#include "trees/FunctionTreeVector.h"

using namespace Eigen;

namespace mrcpp {

/**
 * @brief Lightweight, trivially-copiable metadata for a multicomponent function.
 *
 * This POD accompanies the trees comprising a @ref CompFunction. It holds flags about
 * real/complex storage, conjugation, component counts, user-defined labels, and on-disk
 * layout hints. Arrays have fixed size (4) to simplify MPI packing and shallow copies.
 *
 * @tparam D Spatial dimension of the function (1–3 supported by MRCPP).
 */
template <int D> struct CompFunctionData {
    /** @name Global function descriptors (user-defined) */
    ///@{
    int Ncomp{0};   ///< Number of components actually defined/allocated (0–4).
    int rank{-1};   ///< Rank (index) inside an external vector or basis set.
    int conj{0};    ///< Soft-conjugate flag for algebra (applied to all components).
    int CompFn1{0}; ///< Free integer tag (user purpose).
    int CompFn2{0}; ///< Free integer tag (user purpose).
    int isreal{0};      ///< 1 if component trees are real-valued (`T=double`).
    int iscomplex{0};   ///< 1 if component trees are complex-valued (`T=ComplexDouble`).
    double CompFd1{0.0};///< Free double tag (user purpose).
    double CompFd2{0.0};///< Free double tag (user purpose).
    double CompFd3{0.0};///< Free double tag (user purpose).
    ///@}

    /** @name Per-component user metadata (fixed-size slots 0..3) */
    ///@{
    int n1[4]{0, 0, 0, 0}; ///< Integer label; unequal labels are treated orthogonal in some workflows.
    int n2[4]{0, 0, 0, 0}; ///< Additional integer label (user purpose).
    int n3[4]{0, 0, 0, 0}; ///< Additional integer label (user purpose).
    int n4[4]{0, 0, 0, 0}; ///< Additional integer label (user purpose).

    /**
     * @brief Per-component multiplicative factor.
     *
     * Often used to carry factors like *i* for momentum-like operators without
     * explicitly modifying stored coefficients.
     */
    ComplexDouble c1[4]{{1.0, 0.0}, {1.0, 0.0}, {1.0, 0.0}, {1.0, 0.0}};

    double d1[4]{0.0, 0.0, 0.0, 0.0}; ///< Free double tag (user purpose) per component.
    double d2[4]{0.0, 0.0, 0.0, 0.0}; ///< Free double tag (user purpose) per component.
    double d3[4]{0.0, 0.0, 0.0, 0.0}; ///< Free double tag (user purpose) per component.
    ///@}

    /** @name On-disk/storage layout hints (optional) */
    ///@{
    int type{0};   ///< Serialization type code.
    int order{1};  ///< Polynomial order or filter order hint.
    int scale{0};  ///< Root scale / global scale offset.
    int depth{0};  ///< Max depth.
    int boxes[3] = {0, 0, 0};  ///< Root box tiling (D components used).
    int corner[3] = {0, 0, 0}; ///< Root spatial corner (D components used).
    ///@}

    /** @name Internal runtime fields */
    ///@{
    int shared{0};              ///< 1 if this function uses shared-memory trees.
    int Nchunks[4]{0, 0, 0, 0}; ///< Chunk count for each component (used for MPI shipping).
    ///@}
};

/**
 * @brief Owning pointer wrapper for up to four component trees (real and/or complex).
 *
 * Optionally allocates per-communicator shared memory windows when constructed
 * with @p share = true and MPI shared memory is available (see @ref mpi::comm_share
 * and @ref mpi::shared_memory_size).
 *
 * @tparam D Spatial dimension (1–3).
 */
template <int D> class TreePtr final {
public:
    /**
     * @brief Construct an empty handle.
     * @param share If true and MPI is enabled, create shared-memory windows
     *              for real and complex storage sized per @ref mpi::shared_memory_size (MB).
     */
    explicit TreePtr(bool share)
            : shared_mem_real(nullptr)
            , shared_mem_cplx(nullptr) {
        for (int i = 0; i < 4; i++) real[i] = nullptr;
        for (int i = 0; i < 4; i++) cplx[i] = nullptr;
        is_shared = share;
        if (is_shared and mpi::share_size > 1) {
            // Memory size in MB defined in input. Virtual memory, does not cost anything if not used.
#ifdef MRCPP_HAS_MPI
            this->shared_mem_real = new mrcpp::SharedMemory<double>(mpi::comm_share, mpi::shared_memory_size);
            this->shared_mem_cplx = new mrcpp::SharedMemory<ComplexDouble>(mpi::comm_share, mpi::shared_memory_size);
#endif
        }
    }

    /// Destructor: frees shared windows and any allocated trees.
    ~TreePtr() {
        if (this->shared_mem_real != nullptr) delete this->shared_mem_real;
        if (this->shared_mem_cplx != nullptr) delete this->shared_mem_cplx;
        for (int i = 0; i < 4; i++) {
            if (this->real[i] != nullptr) delete this->real[i];
            if (this->cplx[i] != nullptr) delete this->cplx[i];
            this->real[i] = nullptr;
            this->cplx[i] = nullptr;
        }
    }

    /** @name Metadata forwarding (aliases into @ref data) */
    ///@{
    CompFunctionData<D> data;          ///< Attached function metadata.
    int &Ncomp = data.Ncomp;           ///< Number of active components.
    int &rank = data.rank;             ///< External rank/index tag.
    int &conj = data.conj;             ///< Soft conjugation flag.
    int &isreal = data.isreal;         ///< Real storage flag.
    int &iscomplex = data.iscomplex;   ///< Complex storage flag.
    int &share = data.shared;          ///< Shared-memory flag.
    int *Nchunks = data.Nchunks;       ///< Per-component chunk counts.
    ///@}

    /** True if shared-memory windows were requested/allocated. */
    bool is_shared = false;

    friend class CompFunction<D>;

protected:
    /** Component trees (owned). Slots 0..3 are valid when @ref Ncomp > slot. */
    FunctionTree<D, double> *real[4];        ///< Real components.
    FunctionTree<D, ComplexDouble> *cplx[4]; ///< Complex components.

    /** Optional backing shared-memory windows (one per value type). */
    SharedMemory<double> *shared_mem_real;
    SharedMemory<ComplexDouble> *shared_mem_cplx;
};

/**
 * @brief High-level multicomponent function wrapper on MRCPP trees.
 *
 * A @ref CompFunction manages up to four component trees, either real or complex,
 * and exposes utilities such as allocation, projection, algebraic operations,
 * normalization, conjugation, and data shipping.
 *
 * The class shares its internal state through a `std::shared_ptr<TreePtr<D>>`
 * to enable lightweight copies and move semantics, while retaining clear
 * ownership of the underlying trees.
 *
 * @tparam D Spatial dimension (1–3).
 */
template <int D> class CompFunction {
public:
    /**
     * @name Construction
     * Constructors optionally attach an @ref MultiResolutionAnalysis context,
     * choose component count, and enable shared memory.
     */
    ///@{
    /** @brief Construct empty function bound to @p mra (no components allocated). */
    CompFunction(MultiResolutionAnalysis<D> &mra);
    /** @brief Construct unbound/empty function (MRA set later via allocation). */
    CompFunction();
    /** @brief Construct with @p n1 components (0..4). */
    CompFunction(int n1);
    /**
     * @brief Construct with @p n1 components and shared-memory preference.
     * @param n1    Number of components to allocate (0..4).
     * @param share If true, try to use MPI shared memory for tree storage.
     */
    CompFunction(int n1, bool share);
    /**
     * @brief Construct from metadata @p indata.
     * @param indata Initial metadata (copied).
     * @param alloc  If true, allocate trees according to @p indata.Ncomp.
     */
    CompFunction(const CompFunctionData<D> &indata, bool alloc = false);
    /** @brief Copy constructor: shares underlying pointer (trees may be deep-copied by helpers). */
    CompFunction(const CompFunction<D> &compfunc);
    /** @brief Move constructor. */
    CompFunction(CompFunction<D> &&compfunc);
    /** @brief Copy assignment. */
    CompFunction<D> &operator=(const CompFunction<D> &compfunc);
    ///@}

    /** Virtual destructor. Trees are owned by the shared @ref TreePtr and freed accordingly. */
    virtual ~CompFunction() = default;

    /** @name Raw component access (compatibility aliases) */
    ///@{
    /**
     * @brief Pointer-to-array of real component trees (alias of internal storage).
     * @warning Valid only when @ref isreal() is true.
     */
    FunctionTree<D, double> **CompD;
    /**
     * @brief Pointer-to-array of complex component trees (alias of internal storage).
     * @warning Valid only when @ref iscomplex() is true.
     */
    FunctionTree<D, ComplexDouble> **CompC;
    ///@}

    /** Optional human-readable name. */
    std::string name;

    /** @name Metadata accessors */
    ///@{
    /** @brief Return a copy of the current metadata. */
    CompFunctionData<D> data() const { return func_ptr->data; }
    int Ncomp() const { return func_ptr->data.Ncomp; }         ///< Number of components.
    int rank() const { return func_ptr->data.rank; }           ///< External index/rank.
    int conj() const { return func_ptr->data.conj; }           ///< Soft conjugation flag.
    int isreal() const { return func_ptr->data.isreal; }       ///< Real storage flag.
    int iscomplex() const { return func_ptr->data.iscomplex; } ///< Complex storage flag.
    ///@}

    /** @name Mutators for storage type flags */
    ///@{
    /** @brief Declare that this function stores real-valued components. */
    void defreal() { func_ptr->data.isreal = 1; }
    /** @brief Declare that this function stores complex-valued components. */
    void defcomplex() { func_ptr->data.iscomplex = 1; }
    ///@}

    /** @return 1 if using shared-memory storage. */
    int share() const { return func_ptr->data.shared; }

    /** @return Per-component chunk counts (used for MPI shipping). */
    int *Nchunks() const { return func_ptr->data.Nchunks; }

    /**
     * @brief Copy metadata and optionally allocate components (without copying tree data).
     * @param alloc If true, allocate tree containers for the copied component count.
     * @return A new @ref CompFunction sharing no nodes/coefficients with the source.
     */
    CompFunction paramCopy(bool alloc = false) const;

    /**
     * @brief Integrate the function over the domain.
     * @return Complex integral (real-only functions return real part in `.real()`).
     */
    ComplexDouble integrate() const;

    /**
     * @brief L2 norm of the function.
     * @return \f$\|f\|_2\f$ as a double.
     */
    double norm() const;

    /**
     * @brief Square L2 norm of the function.
     * @return \f$\|f\|_2^2\f$ as a double.
     */
    double getSquareNorm() const;

    /**
     * @brief Allocate @p nalloc component trees.
     * @param nalloc Number of components (0..4). Existing components preserved if possible.
     * @param zero   If true, initialize coefficients to zero.
     */
    void alloc(int nalloc = 1, bool zero = true);

    /**
     * @brief Allocate a single component tree.
     * @param i Component index (0..3).
     */
    void alloc_comp(int i = 0);

    /**
     * @brief Attach an externally created real tree as component @p i.
     * @param tree Ownership is transferred to this object.
     * @param i    Component index (0..3).
     */
    void setReal(FunctionTree<D, double> *tree, int i = 0);

    /**
     * @brief Attach an externally created complex tree as component @p i.
     * @copydetails setReal
     */
    void setCplx(FunctionTree<D, ComplexDouble> *tree, int i = 0);

    /** @brief Set/get external rank/index label. */
    void setRank(int i) { func_ptr->rank = i; };
    const int getRank() const { return func_ptr->rank; };

    /**
     * @brief In-place linear update: @f$f \gets f + c \, g@f$.
     * @param c   Complex scalar.
     * @param inp Addend function (components must be layout-compatible).
     */
    void add(ComplexDouble c, CompFunction<D> inp);

    /**
     * @brief Remove coefficients/nodes below precision @p prec.
     * @param prec Relative (or absolute) precision threshold.
     * @return Number of removed nodes or a non-negative status.
     */
    int crop(double prec);

    /**
     * @brief Multiply the entire function by a complex scalar in-place.
     * @param c Complex factor.
     */
    void rescale(ComplexDouble c);

    /**
     * @brief Release all component trees and reset to empty.
     *
     * Metadata is preserved unless tied to tree content.
     */
    void free();

    /** @return Total memory footprint of nodes (bytes or implementation-defined units). */
    int getSizeNodes() const;

    /** @return Total number of nodes across all component trees. */
    int getNNodes() const;

    /** @brief Flush cached MRA-level data (filters, norms) from component trees. */
    void flushMRAData();

    /** @brief Flush cached function-level data (aux norms, temporaries). */
    void flushFuncData();

    /** @brief Snapshot of the current function metadata (same as @ref data()). */
    CompFunctionData<D> getFuncData() const;

    /** @name Component accessors (non-const/const). */
    ///@{
    FunctionTree<D, double> &real(int i = 0);
    FunctionTree<D, ComplexDouble> &complex(int i = 0);
    const FunctionTree<D, double> &real(int i = 0) const;
    const FunctionTree<D, ComplexDouble> &complex(int i = 0) const;
    ///@}

    /** @name Backwards-compatibility helpers (legacy ComplexFunction interface) */
    ///@{
    void free(int type) { free(); }                ///< Ignored @p type; frees all.
    bool hasReal() const { return isreal(); }      ///< True if real storage is active.
    bool hasImag() const { return iscomplex(); }   ///< True if complex storage is active.
    bool isShared() const { return share(); }      ///< True if shared-memory is active.
    bool conjugate() const { return conj(); }      ///< True if conjugation is requested.
    /** @brief Apply Hermitian adjoint (conjugation + operator-specific flips as implemented). */
    void dagger();
    /** @brief Imaginary component accessor (legacy; identical to @ref real()). */
    FunctionTree<D, double> &imag(int i = 0);
    /** @brief Const imaginary component accessor (legacy; identical to @ref real()). */
    const FunctionTree<D, double> &imag(int i = 0) const;
    ///@}

    /** @brief Shared state (trees + metadata). */
    std::shared_ptr<mrcpp::TreePtr<D>> func_ptr;
};

/** @name Helpers: copying and algebra on @ref CompFunction
 *  Functions operate componentwise on underlying trees and obey precision controls.
 */
///@{
/**
 * @brief Ensure @p out is complex-valued, copying/embedding a real @p inp if needed.
 * @tparam D Dimension.
 */
template <int D> void CopyToComplex(CompFunction<D> &out, const CompFunction<D> &inp);

/** @brief Deep-copy @p inp into *@p out (allocate if needed). */
template <int D> void deep_copy(CompFunction<D> *out, const CompFunction<D> &inp);
/** @brief Deep-copy @p inp into @p out (allocate if needed). */
template <int D> void deep_copy(CompFunction<D> &out, const CompFunction<D> &inp);

/**
 * @brief Compute @f$out = a \, inp\_a + b \, inp\_b@f$ with adaptive precision.
 * @param prec Target precision controlling refinement/cropping.
 * @param conjugate If true, apply soft conjugation to inputs as required.
 */
template <int D>
void add(CompFunction<D> &out, ComplexDouble a, CompFunction<D> inp_a,
         ComplexDouble b, CompFunction<D> inp_b, double prec, bool conjugate = false);

/**
 * @brief Linear combination of many inputs: @f$out = \sum_k c_k \, inp_k@f$.
 * @param c    Coefficients (size must match @p inp).
 * @param inp  Input functions (modified only for temporary workspace).
 * @param prec Target precision.
 * @param conjugate Whether to conjugate inputs (soft).
 */
template <int D>
void linear_combination(CompFunction<D> &out, const std::vector<ComplexDouble> &c,
                        std::vector<CompFunction<D>> &inp, double prec, bool conjugate = false);

/**
 * @brief Pointwise product: @f$out = inp\_a \cdot inp\_b@f$ with refinement control.
 * @param prec      Target precision (relative by default).
 * @param absPrec   If true, treat @p prec as absolute precision.
 * @param useMaxNorms If true, use max norms in error control heuristics.
 * @param conjugate If true, apply soft conjugation to first factor.
 */
template <int D>
void multiply(CompFunction<D> &out, CompFunction<D> inp_a, CompFunction<D> inp_b,
              double prec, bool absPrec = false, bool useMaxNorms = false, bool conjugate = false);

/**
 * @brief Scaled product with loop control: @f$out = coef \cdot inp\_a \cdot inp\_b@f$.
 * @param maxIter Limit iterative refinement steps (-1 for default).
 * @copydetails multiply(CompFunction<D>&,CompFunction<D>,CompFunction<D>,double,bool,bool,bool)
 */
template <int D>
void multiply(double prec, CompFunction<D> &out, double coef,
              CompFunction<D> inp_a, CompFunction<D> inp_b, int maxIter = -1,
              bool absPrec = false, bool useMaxNorms = false, bool conjugate = false);

/** @brief Density from a (possibly complex) function: @f$out = |inp|^2@f$. */
template <int D> void make_density(CompFunction<D> &out, CompFunction<D> inp, double prec);

/** @overload */
template <int D>
void multiply(CompFunction<D> &out, CompFunction<D> inp_a, CompFunction<D> inp_b,
              bool absPrec = false, bool useMaxNorms = false, bool conjugate = false);

/** @brief Multiply by an analytic representable real function @p f. */
template <int D>
void multiply(CompFunction<D> &out, CompFunction<D> &inp_a, RepresentableFunction<D, double> &f,
              double prec, int nrefine = 0, bool conjugate = false);

/** @brief Multiply by an analytic representable complex function @p f. */
template <int D>
void multiply(CompFunction<D> &out, CompFunction<D> &inp_a, RepresentableFunction<D, ComplexDouble> &f,
              double prec, int nrefine = 0, bool conjugate = false);

/** @brief Multiply a single tree by a representable real function. */
template <int D>
void multiply(CompFunction<D> &out, FunctionTree<D, double> &inp_a, RepresentableFunction<D, double> &f,
              double prec, int nrefine = 0, bool conjugate = false);

/** @brief Multiply a single tree by a representable complex function. */
template <int D>
void multiply(CompFunction<D> &out, FunctionTree<D, ComplexDouble> &inp_a, RepresentableFunction<D, ComplexDouble> &f,
              double prec, int nrefine = 0, bool conjugate = false);

/**
 * @brief Complex inner product @f$\langle bra \,|\, ket \rangle@f$.
 * @return Complex inner product consistent with MRCPP normalization.
 */
template <int D> ComplexDouble dot(CompFunction<D> bra, CompFunction<D> ket);

/**
 * @brief Node-wise norm dot helper (diagnostics / preconditioners).
 * @return Real value summarizing node-level contributions.
 */
template <int D> double node_norm_dot(CompFunction<D> bra, CompFunction<D> ket);
///@}

/** @name Projection helpers (3D overloads and templated D)
 *  Project analytic functions onto the multiresolution basis.
 */
///@{
/**
 * @brief Project a real-valued lambda/function @p f onto @p out.
 * @param prec Target precision.
 */
void project(CompFunction<3> &out, std::function<double(const Coord<3> &r)> f, double prec);
/** @brief Real-valued projection (explicit name to avoid overload ambiguities on some compilers). */
void project_real(CompFunction<3> &out, std::function<double(const Coord<3> &r)> f, double prec);
/** @brief Project a complex-valued lambda/function @p f onto @p out. */
void project(CompFunction<3> &out, std::function<ComplexDouble(const Coord<3> &r)> f, double prec);
/** @brief Complex-valued projection (explicit name to avoid overload ambiguities). */
void project_cplx(CompFunction<3> &out, std::function<ComplexDouble(const Coord<3> &r)> f, double prec);

/** @brief Project a representable real function onto @p out. */
template <int D> void project(CompFunction<D> &out, RepresentableFunction<D, double> &f, double prec);
/** @brief Project a representable complex function onto @p out. */
template <int D> void project(CompFunction<D> &out, RepresentableFunction<D, ComplexDouble> &f, double prec);
///@}

/**
 * @brief Orthogonalize @p Ket against @p Bra to precision @p prec (Gram–Schmidt-like).
 * @param prec Target precision controlling projection refinement and cropping.
 */
template <int D> void orthogonalize(double prec, CompFunction<D> &Bra, CompFunction<D> &Ket);

/**
 * @brief Convenience container for 3D composite functions with shared MRA.
 *
 * Provides utilities for distribution and linear-algebra operations across
 * the vector (e.g., rotations and overlap matrices).
 */
class CompFunctionVector : public std::vector<CompFunction<3>> {
public:
    /** @brief Construct a vector with @p N default-initialized functions. */
    CompFunctionVector(int N = 0);

    /** @brief Common MRA pointer for all entries (optional but recommended). */
    MultiResolutionAnalysis<3> *vecMRA;

    /**
     * @brief Distribute internal storage across MPI workers (when enabled).
     *
     * Typically assigns component ownership/ranks and updates metadata so that
     * subsequent parallel operations (add/multiply/dot) can proceed efficiently.
     */
    void distribute();
};

/** @name Vector-level linear algebra and IO utilities */
///@{
/**
 * @brief Apply a unitary (or general) complex rotation @p U in-place: @f$\Phi \gets \Phi U@f$.
 * @param prec Optional precision for intermediate truncations (-1 to keep current).
 */
void rotate(CompFunctionVector &Phi, const ComplexMatrix &U, double prec = -1.0);

/**
 * @brief Apply a rotation @p U producing @p Psi: @f$\Psi \gets \Phi U@f$.
 * @param prec Optional precision for intermediate truncations (-1 to keep current).
 */
void rotate(CompFunctionVector &Phi, const ComplexMatrix &U, CompFunctionVector &Psi, double prec = -1.0);

/**
 * @brief Store per-node coefficient blocks of @p Phi into @p account.
 * @param refTree Reference tree defining the union grid/blocking.
 * @param sizes   Optional fixed block size; -1 to auto-size.
 */
void save_nodes(CompFunctionVector &Phi, mrcpp::FunctionTree<3, double> &refTree,
                BankAccount &account, int sizes = -1);

/**
 * @brief Multiply a vector of functions by a representable function @p f.
 * @param prec    Target precision (-1 to inherit default).
 * @param Func    Optional workspace function (reused).
 * @param nrefine Number of refinement passes (>=0).
 * @param all     If true, apply to all components; else honor component flags.
 * @return Result vector (same size as @p Phi).
 */
CompFunctionVector multiply(CompFunctionVector &Phi, RepresentableFunction<3> &f,
                            double prec = -1.0, CompFunction<3> *Func = nullptr,
                            int nrefine = 1, bool all = false);

/**
 * @brief Set a library-global default MRA used by convenience constructors.
 * @param MRA Non-owning pointer (caller keeps it alive).
 */
void SetdefaultMRA(MultiResolutionAnalysis<3> *MRA);

/**
 * @brief Vectorized inner products: returns @f$\langle Bra_i \,|\, Ket_i \rangle@f$ for all i.
 */
ComplexVector dot(CompFunctionVector &Bra, CompFunctionVector &Ket);

/** @brief Compute the (symmetric) Löwdin overlap matrix @f$S@f$ for @p Phi. */
ComplexMatrix calc_lowdin_matrix(CompFunctionVector &Phi);

/** @brief Overlap matrix of a single set against itself. */
ComplexMatrix calc_overlap_matrix(CompFunctionVector &BraKet);

/** @brief Overlap matrix between two sets @p Bra and @p Ket. */
ComplexMatrix calc_overlap_matrix(CompFunctionVector &Bra, CompFunctionVector &Ket);

/**
 * @brief Pairwise orthogonalization of @p Ket against @p Bra to precision @p prec.
 * @param prec Precision target (relative unless the implementation states otherwise).
 */
void orthogonalize(double prec, CompFunctionVector &Bra, CompFunctionVector &Ket);
///@}

} // namespace mrcpp