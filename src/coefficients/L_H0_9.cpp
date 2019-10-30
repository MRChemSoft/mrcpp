/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2019 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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

// THIS FILE HAS BEEN AUTOGENERATED, DO NOT TOUCH

#include "FilterData.h"

#include <Eigen/Core>


namespace mrcpp {
namespace detail {
auto get_L_H0_9() noexcept -> Eigen::Matrix<double, 10, 10> {
  return (Eigen::Matrix<double, 10, 10>() << 0.707106781186548, 3.81668334298582e-79, 2.06100900521234e-77, 7.63336668597164e-79, 3.05334667438866e-78, 0, -1.10683816946589e-77, 0, 3.5800489757207e-76, -3.81668334298582e-79, -0.612372435695794, 0.353553390593274, -1.85347684843749e-77, 1.53144419137306e-77, -2.0991758386422e-78, 3.76658937410913e-77, 9.27931137763427e-78, -1.33405009972801e-76, -3.07863220153594e-76, 4.08650496463175e-76, -4.1983516772844e-78, -0.684653196881458, 0.176776695296637, -2.59534467323036e-77, 1.06867133603603e-77, -7.09903101795362e-77, -2.13734267207206e-77, 2.55336115645751e-76, 1.46178972036357e-76, -7.90244286165214e-76, 0.233853586673371, 0.405046293650491, -0.522912516583797, 0.0883883476483184, -2.99848185133323e-77, 5.00701146057952e-77, 5.9802657130409e-77, -1.64999991771456e-76, -3.24716262539965e-76, 5.23398484813271e-76, 1.22133866975546e-77, 0.153093108923949, 0.592927061281571, -0.350780380010057, 0.0441941738241592, -9.92337669176313e-78, -8.09136868712993e-77, 1.03050450260617e-77, 5.09908894622905e-76, -2.91976275738415e-77, -0.146575492494482, -0.253876200144874, -0.163876382526586, 0.581703452155821, -0.219863238741723, 0.0220970869120796, 4.68497880351509e-77, -9.4462912738899e-78, -3.21984948522641e-76, 8.22972345831317e-77, 3.12968034124837e-77, -0.0689981317681863, -0.26722861525761, -0.421585548851001, 0.478033079399324, -0.132121363478811, 0.0110485434560398, 1.27477223655726e-76, -7.28986518510291e-77, -2.86442084891086e-76, 0.106977062012728, 0.185289706650491, 0.181798066847189, -0.0566069404148025, -0.534885310063639, 0.3548027758708, -0.0771422564770762, 0.0055242717280199, 1.27858891990025e-77, 1.07296510479689e-76, 3.70218284269624e-76, 0.0394511911654769, 0.152793806371937, 0.301313449619951, 0.204994402553155, -0.528802957971469, 0.246372609862836, -0.0441077726196727, 0.00276213586400995, 3.228914108166e-76, -0.0842790976968415, -0.145975679226991, -0.161531821313557, -0.0637090094933597, 0.216717679791878, 0.39931734961447, -0.455808912293701, 0.163205770906627, -0.0248208301311305, 0.00138106793200498).finished();
}
} // namespace mrcpp
} // namespace detail
