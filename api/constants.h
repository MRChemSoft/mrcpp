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

namespace mrcpp {

const double MachinePrec = 1.0e-15;
const double MachineZero = 1.0e-14;
const int MaxOrder = 41;  ///< Maximum scaling order
const int MaxDepth = 30;  ///< Maximum depth of trees
const int MaxScale = 31;  ///< Maximum scale of trees
const int MinScale = -31; ///< Minimum scale of trees
const int MaxSepRank = 1000;

namespace Axis {
const int None = -1;
const int X = 0;
const int Y = 1;
const int Z = 2;
} // namespace Axis

enum Spin { Paired, Alpha, Beta };
enum FuncType { Legendre, Interpol };
enum SplitType { ExactSplit, NormalSplit, FastSplit };
enum CV_Transform { Forward, Backward };
enum MW_Transform { Compression, Reconstruction };
enum XC_Type { XC_undefined, XC_lda, XC_gga };
enum Traverse { TopDown, BottomUp };
enum Iterator { Lebesgue, Hilbert };

// Math constants
const double pi = 3.1415926535897932384626433832795;
const double root_pi = 1.7724538509055160273;
const double C_x = -0.73855876638202240588; // Dirac exchange constant

} // namespace mrcpp
