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

#include "MultiResolutionAnalysis.h"
#include "core/FilterCache.h"
#include "core/InterpolatingBasis.h"
#include "core/LegendreBasis.h"
#include "utils/Printer.h"
#include "utils/math_utils.h"

namespace mrcpp {

/** @returns New MultiResolutionAnalysis (MRA) object
 *
 *  @brief Contructs a MultiResolutionAnalysis object composed of computational domain (world) and a polynomial basis (Multiwavelets)
 *
 *  @param[in] bb: 2-element integer array [Lower, Upper] defining the bounds for a BoundingBox object representing the computational domain
 *  @param[in] order: (integer) Maximum polynomial order of the multiwavelet basis, 
 *  immediately used in the constructor of an InterPolatingBasis object which becomes an attribute of the MRA
 *  @param[in] maxDepth: (integer) Exponent of the node refinement in base 2, relative to root scale. 
 *  In other words, it is the maximum amount of refinement that we allow in a node, in other to avoid overflow of values.
 *
 *  @details Contructor of the MultiResolutionAnalysis class from scratch, without requiring any pre-existing complex structure. 
 *  The constructor calls the InterpolatingBasis basis constructor to generate the MultiWavelets basis of functions,
 *  then the BoundingBox constructor to create the computational domain. The constructor then checks if the generated node depth, or
 *  node refinement is beyond the root scale or the maximum depth allowed, in which case it will abort the process.
 *  Otherwise, the process goes on to setup the filters with the class' setupFilter method.
 */
template <int D>
MultiResolutionAnalysis<D>::MultiResolutionAnalysis(std::array<int, 2> bb, int order, int depth)
        : maxDepth(depth)
        , basis(InterpolatingBasis(order))
        , world(bb) {
    if (getMaxDepth() > MaxDepth) MSG_ABORT("Beyond MaxDepth");
    if (getMaxScale() > MaxScale) MSG_ABORT("Beyond MaxScale");
    setupFilter();
}

/** @returns New MultiResolutionAnalysis (MRA) object
 *
 *  @brief Contructs a MultiResolutionAnalysis object composed of computational domain (world) and a polynomial basis (Multiwavelets) from a pre-existing BoundingBox object
 *
 *  @param[in] bb: BoundingBox object representing the computational domain
 *  @param[in] order: (integer) Maximum polynomial order of the multiwavelet basis, 
 *  immediately used in the constructor of an InterPolatingBasis object which becomes an attribute of the MRA
 *  @param[in] maxDepth: (integer) Exponent of the node refinement in base 2, relative to root scale. 
 *  In other words, it is the maximum amount of refinement that we allow in a node, in other to avoid overflow of values.
 *
 *  @details Contructor of the MultiResolutionAnalysis class from a BoundingBox object. For more details see the first constructor.
 */
template <int D>
MultiResolutionAnalysis<D>::MultiResolutionAnalysis(const BoundingBox<D> &bb, int order, int depth)
        : maxDepth(depth)
        , basis(InterpolatingBasis(order))
        , world(bb) {
    if (getMaxDepth() > MaxDepth) MSG_ABORT("Beyond MaxDepth");
    if (getMaxScale() > MaxScale) MSG_ABORT("Beyond MaxScale");
    setupFilter();
}

/** @returns New MultiResolutionAnalysis (MRA) object
 *
 *  @brief Copy constructor for a MultiResolutionAnalysis object composed of computational domain (world) and a polynomial basis (Multiwavelets)
 *
 *  @param[in] mra: Pre-existing MRA object
 *
 *  @details Copy a MultiResolutionAnalysis object without modifying the original. For more details see the first constructor.
 */
template <int D>
MultiResolutionAnalysis<D>::MultiResolutionAnalysis(const MultiResolutionAnalysis<D> &mra)
        : maxDepth(mra.maxDepth)
        , basis(mra.basis)
        , world(mra.world) {
    if (getMaxDepth() > MaxDepth) MSG_ABORT("Beyond MaxDepth");
    if (getMaxScale() > MaxScale) MSG_ABORT("Beyond MaxScale");
    setupFilter();
}

/** @returns New MultiResolutionAnalysis object
 *
 *  @brief Constructor for a MultiResolutionAnalysis object from a pre-existing BoundingBox (computational domain) and a ScalingBasis (Multiwavelet basis) objects
 *
 * @param[in] bb: Computational domain as a BoundingBox object, taken by constant reference
 * @param[in] sb: Polynomial basis (MW) as a ScalingBasis object
 * @param[in] depth: Maximum allowed resolution depth, relative to root scale
 *
 *  @details Creates a MRA object from pre-existing BoundingBox and ScalingBasis objects. These objects are taken as reference. For more details about the constructor itself, see the first constructor.
 */
template <int D>
MultiResolutionAnalysis<D>::MultiResolutionAnalysis(const BoundingBox<D> &bb, const ScalingBasis &sb, int depth)
        : maxDepth(depth)
        , basis(sb)
        , world(bb) {
    if (getMaxDepth() > MaxDepth) MSG_ABORT("Beyond MaxDepth");
    if (getMaxScale() > MaxScale) MSG_ABORT("Beyond MaxScale");
    setupFilter();
}

/** @returns Boolean value
 *
 *  @brief Equality operator for the MultiResolutionAnalysis class, returns true if both MRAs have the same polynomial basis, computational domain and maximum depth, and false otherwise
 *
 * @param[in] mra: MRA object, taken by constant reference
 *
 *  @details Equality operator for the MultiResolutionAnalysis class, returns true if both MRAs have the same polynomial basis represented by a BoundingBox object, computational domain (ScalingBasis object) and maximum depth (integer), and false otherwise.
 *  Computations on different MRA cannot be combined, this operator can be used to make sure that the multiple MRAs are compatible.
 *  For more information about the meaning of equality for BoundingBox and ScalingBasis objets, see their respective classes.
 */
template <int D> bool MultiResolutionAnalysis<D>::operator==(const MultiResolutionAnalysis<D> &mra) const {
    if (this->basis != mra.basis) return false;
    if (this->world != mra.world) return false;
    if (this->maxDepth != mra.maxDepth) return false;
    return true;
}

/** @returns Boolean value
 *
 *  @brief Inequality operator for the MultiResolutionAnalysis class, returns false if both MRAs have the same polynomial basis, computational domain and maximum depth, and true otherwise
 *
 * @param[in] mra: MRA object, taken by constant reference
 *
 *  @details Inequality operator for the MultiResolutionAnalysis class, returns true if both MRAs have the same polynomial basis represented by a BoundingBox object, computational domain (ScalingBasis object) and maximum depth (integer), and false otherwise.
 *  Opposite of the == operator.
 *  For more information about the meaning of equality for BoundingBox and ScalingBasis objets, see their respective classes.
 */
template <int D> bool MultiResolutionAnalysis<D>::operator!=(const MultiResolutionAnalysis<D> &mra) const {
    if (this->basis != mra.basis) return true;
    if (this->world != mra.world) return true;
    if (this->maxDepth != mra.maxDepth) return true;
    return false;
}

/** @returns void
 *
 *  @brief Displays the MRA's attributes in the outstream defined in the Printer class
 *
 *  @details This function displays the attributes of the MRA in the using the Printer class. 
 *  By default, the Printer class writes all information in the output file, not the terminal.
 *  
 */
template <int D> void MultiResolutionAnalysis<D>::print() const {
    print::separator(0, ' ');
    print::header(0, "MultiResolution Analysis");
    println(0, this->basis);
    print::separator(0, '-');
    println(0, this->world);
    print::separator(0, '=', 2);
}

/** @returns void
 *
 *  @brief TODO: Sets up the filters
 *
 *  @details TBD
 *  
 */
template <int D> void MultiResolutionAnalysis<D>::setupFilter() {
    getLegendreFilterCache(lfilters);
    getInterpolatingFilterCache(ifilters);
    int k = this->basis.getScalingOrder();
    int type = this->basis.getScalingType();
    switch (type) {
        case Legendre:
            this->filter = &lfilters.get(k);
            break;
        case Interpol:
            this->filter = &ifilters.get(k);
            break;
        default:
            MSG_ERROR("Invalid scaling basis selected.")
    }
}

/** @returns double (double-precision floating point number) 
 *
 *  @brief Computes the difference between the lower and upper bounds of the computational domain 
 *
 *  @details This function displays the attributes of the MRA in the using the Printer class. 
 *  By default, the Printer class writes all information in the output file, not the terminal.
 *  
 */
template <int D> double MultiResolutionAnalysis<D>::calcMaxDistance() const {
    const Coord<D> &lb = getWorldBox().getLowerBounds();
    const Coord<D> &ub = getWorldBox().getUpperBounds();
    return math_utils::calc_distance<D>(lb, ub);
}

template class MultiResolutionAnalysis<1>;
template class MultiResolutionAnalysis<2>;
template class MultiResolutionAnalysis<3>;

} // namespace mrcpp
