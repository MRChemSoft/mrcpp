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

#include <array>
#include <iomanip>
#include <ostream>

#include "NodeIndex.h"
#include "utils/details.h"

#include "MRCPP/mrcpp_declarations.h"

namespace mrcpp {

/** 
 * @class BoundingBox
 * @tparam D Spatial dimension (1, 2, or 3)
 *
 * @brief Class defining the computational domain
 *
 * @details The computational domain is made up of a collection of D-dimensional
 * boxes on a particular length scale \f$ n \f$. The size of each box is then
 * \f$ [2^{-n}]^D \f$, i.e. higher scale means smaller boxes, and the scale
 * may be negative. The number of boxes can be different in each dimension
 * \f$ [n_x, n_y, \dots] \f$, but they must all be on the same scale (size).
 * Box translations relative to the world origin _must_ be an integer
 * multiple of the given scale size \f$ 2^{-n} \f$.
 */
template <int D> class BoundingBox {
public:
    /** 
     * @brief Constructor for BoundingBox object
     * @param box [lower, upper] bound in all dimensions
     *
     * @details Creates a box with appropriate root scale and scaling
     * factor to fit the given bounds, which applies to _all_ dimensions.
     * Root scale is chosen such that the scaling factor becomes 1 < sfac < 2.
     *
     * @note Limitations: Box must be _either_ [0,L] _or_ [-L,L], with L a positive integer.
     * This is the most general constructor, which will create a world with no periodic boundary conditions.
     */
    explicit BoundingBox(std::array<int, 2> box);

    /** 
     * @brief Constructor for BoundingBox object
     * @param n Length scale, default 0
     * @param l Corner translation, default [0, 0, ...]
     * @param nb Number of boxes, default [1, 1, ...]
     * @param sf Scaling factor, default [1.0, 1.0, ...]
     * @param pbc Periodic boundary conditions, default false
     *
     * @details Creates a box with given parameters. The parameter n defines the length scale, which, together with sf, determines the unit length of each side of the boxes by \f$ [2^{-n}]^D \f$.
     * The parameter l defines the corner translation of the lower corner of the box relative to the world origin.
     * The parameter nb defines the number of boxes in each dimension.
     * The parameter sf defines the scaling factor, which determines the box translations around the origin, i.e. the amount of boxes around origin.
     * The parameter pbc defines whether the world is periodic or not. In this constructor this value makes all dimensions periodic.
     * This constructor is used for work in periodic systems.
     */
    explicit BoundingBox(int n = 0,
                         const std::array<int, D> &l = {},
                         const std::array<int, D> &nb = {},
                         const std::array<double, D> &sf = {},
                         bool pbc = false);

    /** 
     * @brief Constructor for BoundingBox objec
     * @param idx index of the lower corner of the box
     * @param nb Number of boxes, default [1, 1, ...]
     * @param sf Scaling factor, default [1.0, 1.0, ...]
     *
     * @details Creates a box with given parameters
     * The parameter idx defines the index of the lower corner of the box relative to the world origin.
     * The parameter nb defines the number of boxes in each dimension.
     * The parameter sf defines the scaling factor, which determines the box translations around the origin, i.e. the amount of boxes around origin.
     * This constructor creates a world with no periodic boundary conditions.
     */
    explicit BoundingBox(const NodeIndex<D> &idx,
                         const std::array<int, D> &nb = {},
                         const std::array<double, D> &sf = {});

    /** 
     * @brief Constructor for BoundingBox object
     * @param sf Scaling factor, default [1.0, 1.0, ...]
     * @param pbc Periodic boundary conditions, default true
     *
     * @details Creates a box with given parameters.
     * The parameter sf defines the scaling factor, which determines the box translations around the origin, i.e. the amount of boxes around origin.
     * The parameter pbc defines whether the world is periodic or not. In this constructor this value makes all dimensions periodic.
     * This construtor is used for work in periodic systems.
     */
    explicit BoundingBox(const std::array<double, D> &sf, bool pbc = true);

    /** 
     * @brief Constructor for BoundingBox object
     * @param sf Scaling factor, default [1.0, 1.0, ...]
     * @param pbc Periodic boundary conditions, default true
     *
     * @details Creates a box with given parameters.
     * The parameter sf defines the scaling factor, which determines the box translations around the origin, i.e. the amount of boxes around origin.
     * The parameter pbc defines whether the world is periodic or not. In this constructor this value makes specific dimensions periodic.
     * This is used for work in periodic systems.
     */
    BoundingBox(const std::array<double, D> &sf, std::array<bool, D> pbc);

    /** 
     * @brief Copy constructor for BoundingBox object
     * @param box Other BoundingBox object
     *
     * @details Creates a box identical to the input box paramter.
     * This constructor uses all parameters from the other BoundingBox to create a new one.
     */
    BoundingBox(const BoundingBox<D> &box);

    /** 
     * @brief Assignment operator overload for BoundingBox object
     * @param box Other BoundingBox object
     * @details Allocates all parameters in this BoundingBox to be that of the other BoundingBox.
     * @return New BoundingBox object
     */
    BoundingBox<D> &operator=(const BoundingBox<D> &box);

    virtual ~BoundingBox() = default;

    /**
     * @brief Equality: same corner index and per-dimension box counts
     * @param box Other BoundingBox object
     * @return True if equal, false otherwise
     */
    inline bool operator==(const BoundingBox<D> &box) const;
    /**
     * @brief Inequality: differs in corner index or in any per-dimension box count
     * @param box Other BoundingBox object
     * @return True if not equal, false otherwise
     */
    inline bool operator!=(const BoundingBox<D> &box) const;

    /** 
     * @brief Fetches a NodeIndex object from a given box index
     * @param bIdx The index of the box we want to fetch the cell index from
     *
     * @details During the adaptive refinement, each original box will contain an increasing number of smaller cells,
     * each of which will be part of a specific node in the tree. These cells are divided adaptivelly. This function returns the NodeIndex
     * object of the cell at the lower back corner of the box object indexed by bIdx.
     * 
     * @return The NodeIndex object of the index given as it is in the Multiresolutoin analysis
     * @note Specialized for D=1 below
     */
    NodeIndex<D> getNodeIndex(int bIdx) const;

    /** 
     * @brief Fetches the index of a box from a given coordinate
     * @param r D-dimensional array representaing a coordinate in the simulation box
     * @return The index value of the boxes in the position given as it is in the generated world
     * @note Specialized for D=1 below
     */
    int getBoxIndex(Coord<D> r) const;

    /** 
     * @brief Fetches the index of a box from a given NodeIndex
     * @param nIdx NodeIndex object, representing the node and its index in the adaptive tree
     *
     * @details During the multiresolution analysis the boxes will be divided into smaller boxes, which means that each individual box will be part of a specific node in the tree.
     * Each node will get its own index value, but will still be part of one of the original boxes of the world.
     * 
     * @return The index value of the boxes in which the NodeIndex object is mapping to
     * @note Specialized for D=1 below
     */
    int getBoxIndex(NodeIndex<D> nIdx) const;

    int size() const { return this->totBoxes; }                     ///< @return Total number of boxes
    /**
     * @param d Dimension index
     * @return Number of boxes along dimension @p d
     */
    int size(int d) const { return this->nBoxes[d]; }
    int getScale() const { return this->cornerIndex.getScale(); }   ///< @return Root scale \(n\)

    /**
     * @param d Dimension index
     * @return Scaling factor to scale this box by along dimension @p d
     */
    double getScalingFactor(int d) const { return this->scalingFactor[d]; }
    /**
     * @param d Dimension index
     * @return Unit length along dimension @p d
     */
    double getUnitLength(int d) const { return this->unitLengths[d]; }
    /**
     * @param d Dimension index
     * @return Box length along dimension @p d
     */
    double getBoxLength(int d) const { return this->boxLengths[d]; }
    /**
     * @param d Dimension index
     * @return Lower bound of this box coordinates along dimension @p d
     */
    double getLowerBound(int d) const { return this->lowerBounds[d]; }
    /**
     * @param d Dimension index
     * @return Upper bound of this box coordinates along dimension @p d
     */
    double getUpperBound(int d) const { return this->upperBounds[d]; }

    bool isPeriodic() const { return details::are_any(this->periodic, true); } ///< @return Is any dimension periodic?
    const std::array<bool, D> &getPeriodic() const { return this->periodic; }  ///< @return Periodicity flags per dimension

    const Coord<D> &getUnitLengths() const { return this->unitLengths; }                    ///< @return The unit lengths
    const Coord<D> &getBoxLengths() const { return this->boxLengths; }                      ///< @return The box lengths
    const Coord<D> &getLowerBounds() const { return this->lowerBounds; }                    ///< @return The lower bounds of the coordinates of this box
    const Coord<D> &getUpperBounds() const { return this->upperBounds; }                    ///< @return The upper bounds of the coordinates of this box
    const NodeIndex<D> &getCornerIndex() const { return this->cornerIndex; }                ///< @return The corner index
    const std::array<double, D> &getScalingFactors() const { return this->scalingFactor; }  ///< @return The scaling factors to scale this box by

    /**
     * @brief Stream output operator
     * @param o Output stream
     * @param box BoundingBox object
     * @return Reference to output stream
     */
    friend std::ostream &operator<<(std::ostream &o, const BoundingBox<D> &box) { return box.print(o); }

protected:
    // Fundamental parameters
    NodeIndex<D> cornerIndex;               ///< Index defining the lower corner of the box
    std::array<int, D> nBoxes{};            ///< Number of boxes in each dim, last entry total
    std::array<double, D> scalingFactor{};  ///< Scaling factors to scale this box by, per dimension
    std::array<bool, D> periodic{};         ///< Sets which dimension has Periodic boundary conditions

    // Derived parameters

    int totBoxes{1};       ///< Number of total boxes
    Coord<D> unitLengths;  ///< 1/2^initialScale
    Coord<D> boxLengths;   ///< Total length (unitLength times nBoxes)
    Coord<D> lowerBounds;  ///< Box lower bound (not real)
    Coord<D> upperBounds;  ///< Box upper bound (not real)

    /** 
     * @brief Sets the number of boxes in each dimension
     * @param nb Number of boxes, default [1, 1, ...]
     *
     * @details For each dimentions D it sets the number of boxes in that dimension in the nBoxes array and the total amount of boxes in the world in the totBoxes variable.
     * This just sets counters for the number of boxes in each dimension.
     */
    void setNBoxes(const std::array<int, D> &nb = {});

    /** 
     * @brief Computes and sets all derived parameters
     * 
     * @details For all parameters that have been initialized in the constructor,
     * this function will compute the necessary derived parameters in each dimension.
     * The unit length is set to \a sfac \f$ \cdot 2^{-n} \f$ where \a sfac is the scaling factor (default 1.0) and n is the length scale.
     * The unit length is the base unit which is used for the size and positioning of the boxes around origin.
     * The boxLength is the total length of the box in each dimension, which is the unit length times the number of boxes in that dimension.
     * The lowerBound is computed from the index of the lower corner of the box and the unit length.
     * The upperBound is computed to be the lower corner plus the total length in that dimension.
     */
    void setDerivedParameters();

    /** 
     * @brief Sets the scaling factors in each dimension
     * @param sf Scaling factor, default [1.0, 1.0, ...]
     * 
     * @details This checks that the sf variable has sane values before assigning it to the member variable scalingFactor.
     */
    void setScalingFactors(const std::array<double, D> &sf);

    /** 
     * @brief Sets whether all dimensions are periodic
     * @param pbc Boolean which is used to set all dimension to either periodic or not
     *
     * @details This fills in the periodic array with the values from the input.
     */
    void setPeriodic(std::array<bool, D> periodic);

    /** 
     * @brief Sets which dimensions are periodic
     * @param pbc D-dimensional array holding boolean values for each dimension
     *
     * @details This fills in the periodic array with the values from the input array.
     */
    void setPeriodic(bool periodic);

    /**
     * @brief Prints information about the BoundinBox object
     * @param o Output stream variable which will be used to print the information
     * 
     * @details A function which prints information about the BoundingBox object.
     * 
     * @return The output stream variable
     */
    std::ostream &print(std::ostream &o) const;
};

// Inline comparisons

template <int D> bool BoundingBox<D>::operator==(const BoundingBox<D> &box) const {
    if (getCornerIndex() != box.getCornerIndex()) return false;
    for (int d = 0; d < D; d++) {
        if (this->size(d) != box.size(d)) return false;
    }
    return true;
}

template <int D> bool BoundingBox<D>::operator!=(const BoundingBox<D> &box) const {
    if (getCornerIndex() != box.getCornerIndex()) return true;
    for (int d = 0; d < D; d++) {
        if (this->size(d) != box.size(d)) return true;
    }
    return false;
}

} // namespace mrcpp