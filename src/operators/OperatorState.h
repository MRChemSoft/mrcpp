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

/**
 * @file
 * @brief Lightweight state holder used during operator application.
 *
 * @details
 * The operator application kernels (e.g., convolution/derivative calculators)
 * are performance-critical and multi-threaded. To avoid sharing mutable
 * state between threads, this helper encapsulates:
 *  - pointers to the current *source* (g) and *destination* (f) MW nodes,
 *  - precomputed size/stride quantities (`kp1`, `kp1_d`, …),
 *  - addresses of coefficient blocks for selected components (ft/gt),
 *  - temporary scratch buffers laid out for cache-friendly sweeps,
 *  - and small per-call metadata such as the maximum index offset
 *    (`maxDeltaL`) between the active nodes.
 *
 * It is deliberately simple (POD-like) and header-only to enable aggressive
 * inlining by the compiler.
 */

#pragma once

#include <vector>

#include <Eigen/Core>

#include "trees/MWNode.h"
#include "utils/math_utils.h"

namespace mrcpp {

/**
 * @def GET_OP_IDX(FT, GT, ID)
 * @brief Build a 2-bit operator index for dimension @p ID from component flags.
 *
 * @details
 * Encodes the *from* (FT) and *to* (GT) component bits at position @p ID into
 * an index in the set {0,1,2,3}:
 * \f[
 *   \mathrm{idx} = 2 \cdot \big( (GT \gg ID) \& 1 \big)
 *                +      \big( (FT \gg ID) \& 1 \big).
 * \f]
 * This compact index is used to select per-dimension operator blocks.
 */
#define GET_OP_IDX(FT, GT, ID) (2 * ((GT >> ID) & 1) + ((FT >> ID) & 1))

/**
 * @class OperatorState
 * @brief Thread-local state for applying an MW operator to node data.
 *
 * @tparam D Spatial dimension of the node (1–3).
 * @tparam T Coefficient value type (e.g., double or std::complex<double>).
 *
 * @details
 * The class provides:
 *  - Binding of a *g-node* (source) at construction time.
 *  - Late binding of an *f-node* (destination) and its @ref NodeIndex.
 *  - Selection of component blocks (ft/gt) via @ref setFComponent and
 *    @ref setGComponent, exposing the corresponding coefficient slices.
 *  - Access to alternating scratch buffers arranged as
 *    `aux[0] = f-comp`, `aux[1..D-1]` ping-pong across `scr1`/`scr2`,
 *    and `aux[D] = g-comp`.
 *
 * The scratch layout avoids reallocation and reduces cache conflicts during
 * dimension-by-dimension tensor sweeps.
 */
template <int D, typename T> class OperatorState final {
public:
    /**
     * @brief Construct with a source node and a raw scratch buffer.
     *
     * @param gn   Source (g) node whose coefficients are read.
     * @param scr1 Pointer to a scratch buffer of at least `kp1_d` elements.
     *
     * @details
     * Two scratch regions are interleaved: `scr1` and `scr2 = scr1 + kp1_d`.
     * For each interior dimension `i=1..D-1`, the buffer alternates between
     * these two regions by parity of `i` to enable out-of-place 1D transforms.
     */
    OperatorState(MWNode<D, T> &gn, T *scr1)
            : gNode(&gn) {
        this->kp1      = this->gNode->getKp1();                // basis points per dim
        this->kp1_d    = this->gNode->getKp1_d();              // total points (kp1^D)
        this->kp1_2    = math_utils::ipow(this->kp1, 2);       // kp1^2
        this->kp1_dm1  = math_utils::ipow(this->kp1, D - 1);   // kp1^(D-1)
        this->gData    = this->gNode->getCoefs();
        this->maxDeltaL = -1;

        T *scr2 = scr1 + this->kp1_d;

        // Assign alternating aux buffers for interior dimensions
        for (int i = 1; i < D; i++) {
            if (IS_ODD(i)) {
                this->aux[i] = scr2;
            } else {
                this->aux[i] = scr1;
            }
        }
    }

    /**
     * @brief Convenience ctor: scratch storage provided as a std::vector.
     * @param gn   Source (g) node.
     * @param scr1 Vector whose data pointer is used as scratch.
     *
     * @warning The vector must outlive the OperatorState.
     */
    OperatorState(MWNode<D, T> &gn, std::vector<T> scr1)
            : OperatorState(gn, scr1.data()) {}

    /**
     * @brief Bind the destination (f) node and cache its coefficient pointer.
     */
    void setFNode(MWNode<D, T> &fn) {
        this->fNode = &fn;
        this->fData = this->fNode->getCoefs();
    }

    /**
     * @brief Bind the destination node index and update @ref maxDeltaL.
     * @param idx Destination node index in the tree.
     *
     * @details
     * The maximum level shift \f$\max_d |f_l[d] - g_l[d]|\f$ is used to
     * select scale-dependent operator stencils/bandwidths.
     */
    void setFIndex(NodeIndex<D> &idx) {
        this->fIdx = &idx;
        calcMaxDeltaL();
    }

    /**
     * @brief Select the source (g) component and expose its coefficient slice.
     * @param gt Component bitfield (typically 0/1 per dimension).
     *
     * @details Offsets the base pointer by `gt * kp1_d` and stores it in
     * `aux[D]`, which operator kernels read as the final stage input.
     */
    void setGComponent(int gt) {
        this->aux[D] = this->gData + gt * this->kp1_d;
        this->gt = gt;
    }

    /**
     * @brief Select the destination (f) component and expose its coefficient slice.
     * @param ft Component bitfield (typically 0/1 per dimension).
     *
     * @details Offsets the base pointer by `ft * kp1_d` and stores it in
     * `aux[0]`, which operator kernels use as the first stage buffer.
     */
    void setFComponent(int ft) {
        this->aux[0] = this->fData + ft * this->kp1_d;
        this->ft = ft;
    }

    /**
     * @brief Maximum level difference between the bound f/g nodes.
     */
    int getMaxDeltaL() const { return this->maxDeltaL; }

    /**
     * @brief Build a compact operator index for dimension @p i (0..D-1).
     *
     * @details Uses @ref GET_OP_IDX on the currently bound @ref ft and @ref gt.
     */
    int getOperIndex(int i) const { return GET_OP_IDX(this->ft, this->gt, i); }

    /**
     * @brief Access the array of auxiliary data pointers used by kernels.
     * @return `aux[0] = f-comp`, `aux[1..D-1]` scratch, `aux[D] = g-comp`.
     */
    T **getAuxData() { return this->aux; }

    /**
     * @brief Access per-dimension operator data blocks (set by calculators).
     */
    double **getOperData() { return this->oData; }

    // Calculator kernels are declared as friends to allow fast access.
    friend class ConvolutionCalculator<D, T>;
    friend class DerivativeCalculator<D, T>;

private:
    // Current component selectors (bitfields)
    int ft{0};
    int gt{0};

    // Geometry / thresholds
    int maxDeltaL;     ///< max_d |f_l[d] - g_l[d]|; computed in calcMaxDeltaL()
    double fThreshold; ///< (optional) threshold for f (may be set by calculators)
    double gThreshold; ///< (optional) threshold for g (may be set by calculators)

    // Shorthands derived from the bound node
    int kp1;     ///< #points per dimension
    int kp1_2;   ///< kp1^2
    int kp1_d;   ///< kp1^D (total points in a component block)
    int kp1_dm1; ///< kp1^(D-1)

    // Bound nodes and indices
    MWNode<D, T> *gNode{nullptr};
    MWNode<D, T> *fNode{nullptr};
    NodeIndex<D> *fIdx{nullptr};

    // Data pointers
    T *aux[D + 1]{}; ///< [0]=f-comp, [1..D-1]=scratch, [D]=g-comp
    T *gData{nullptr};
    T *fData{nullptr};
    double *oData[D]{}; ///< Per-dimension operator-specific metadata

    /// @brief Compute @ref maxDeltaL from the currently bound f/g nodes.
    void calcMaxDeltaL() {
        const auto &gl = this->gNode->getNodeIndex();
        const auto &fl = *this->fIdx;
        int max_dl = 0;
        for (int d = 0; d < D; d++) {
            int dl = std::abs(fl[d] - gl[d]);
            if (dl > max_dl) { max_dl = dl; }
        }
        this->maxDeltaL = max_dl;
    }
};

} // namespace mrcpp