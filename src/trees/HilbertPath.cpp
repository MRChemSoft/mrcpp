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

#include "HilbertPath.h"

namespace mrcpp {
// clang-format off
template<>
const short int HilbertPath<1>::pTable[1][8] = {
    {0,0,-1,-1,-1,-1,-1,-1}
};

template<>
const int HilbertPath<1>::zTable[1][8] = {
    {0,1,-1,-1,-1,-1,-1,-1}
};

template<>
const int HilbertPath<1>::hTable[1][8] = {
    {0,1,-1,-1,-1,-1,-1,-1}
};

template<>
const short int HilbertPath<2>::pTable[4][8] = {
    {1,0,0,3,-1,-1,-1,-1},
    {0,1,1,2,-1,-1,-1,-1},
    {3,2,2,1,-1,-1,-1,-1},
    {2,3,3,0,-1,-1,-1,-1}
};

template<>
const int HilbertPath<2>::zTable[4][8] = {
    {0,2,3,1,-1,-1,-1,-1},
    {0,1,3,2,-1,-1,-1,-1},
    {3,1,0,2,-1,-1,-1,-1},
    {3,2,0,1,-1,-1,-1,-1}
};

template<>
const int HilbertPath<2>::hTable[4][8] = {
    {0,3,1,2,-1,-1,-1,-1},
    {0,1,3,2,-1,-1,-1,-1},
    {2,1,3,0,-1,-1,-1,-1},
    {2,3,1,0,-1,-1,-1,-1}
};

template<>
const short int HilbertPath<3>::pTable[12][8] = {
    { 1, 2, 2, 9, 9, 8, 8, 4},
    { 2, 0, 0, 7, 7, 3, 3,11},
    { 0, 1, 1, 5, 5,10,10, 6},
    { 4, 5, 5, 6, 6,11,11, 1},
    { 5, 3, 3,10,10, 0, 0, 8},
    { 3, 4, 4, 2, 2, 7, 7, 9},
    { 7, 8, 8, 3, 3, 2, 2,10},
    { 8, 6, 6, 1, 1, 9, 9, 5},
    { 6, 7, 7,11,11, 4, 4, 0},
    {10,11,11, 0, 0, 5, 5, 7},
    {11, 9, 9, 4, 4, 6, 6, 2},
    { 9,10,10, 8, 8, 1, 1, 3}
};

template<>
const int HilbertPath<3>::zTable[12][8] = {
    {0,2,6,4,5,7,3,1},
    {0,4,5,1,3,7,6,2},
    {0,1,3,2,6,7,5,4},
    {3,1,5,7,6,4,0,2},
    {3,7,6,2,0,4,5,1},
    {3,1,0,2,6,4,5,7},
    {5,7,3,1,0,2,6,4},
    {5,1,0,4,6,2,3,7},
    {5,4,6,7,3,2,0,1},
    {6,4,0,2,3,1,5,7},
    {6,2,3,7,5,1,0,4},
    {6,7,5,4,0,1,3,2}
};

template<>
const int HilbertPath<3>::hTable[12][8] = {
    {0,7,1,6,3,4,2,5},
    {0,3,7,4,1,2,6,5},
    {0,1,3,2,7,6,4,5},
    {6,1,7,0,5,2,4,3},
    {4,7,3,0,5,6,2,1},
    {2,1,3,0,5,6,4,7},
    {4,3,5,2,7,0,6,1},
    {2,1,5,6,3,0,4,7},
    {6,7,5,4,1,0,2,3},
    {2,5,3,4,1,6,0,7},
    {6,5,1,2,7,4,0,3},
    {4,5,7,6,3,2,0,1}
};

template <int D>
short int HilbertPath<D>::getChildPath(int hIdx) const {
    return pTable[this->path][hIdx];
}

template <int D>
int HilbertPath<D>::getZIndex(int hIdx) const {
    return zTable[this->path][hIdx];
}

template <int D>
int HilbertPath<D>::getHIndex(int zIdx) const {
    return hTable[this->path][zIdx];
}

// clang-format on

template class HilbertPath<1>;
template class HilbertPath<2>;
template class HilbertPath<3>;

} // namespace mrcpp
