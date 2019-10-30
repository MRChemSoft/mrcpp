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
auto get_I_G0_15() noexcept -> Eigen::Matrix<double, 16, 16> {
  return (Eigen::Matrix<double, 16, 16>() << 0.188608650494982, -0.289421601111292, 0.155101741373147, 0.113854087356149, -0.281265724910027, 0.169817550450313, 0.128759983378058, -0.265737480503766, 0.0265772582319394, 0.243589202566375, -0.0428623943630423, -0.220089018412713, -0.0875009055097845, 0.0774761922894686, 0.126286060935075, 0.0925910951883618, -0.165099706132016, 0.285770484410546, -0.232967923354118, 0.0267731704956894, 0.202613235556883, -0.276505083934221, 0.104176499571325, 0.17677097342124, -0.242673652508637, -0.0370884565211736, 0.235456542831691, 0.0489297156225168, -0.167671157444739, -0.160873012189309, -0.065655962835585, -0.00861116283026632, -0.139710245186262, 0.266063650630646, -0.279372120134331, 0.160860284930976, 0.0463727901942574, -0.233140157766698, 0.262197434221179, -0.0722265760628254, -0.192985567276492, 0.216997760523271, 0.0846716838309248, -0.2013846867567, -0.128521992292311, 0.0558879427749035, 0.122644192630747, 0.0925090309235043, 0.112782357646123, -0.231520362316204, 0.288197010796437, -0.255690582258724, 0.127835904621231, 0.0630592934213205, -0.231498703858627, 0.255593988189245, -0.067836991349937, -0.193190060894362, 0.189835955273826, 0.12566024819396, -0.136089083257765, -0.16852805127212, -0.0723720215171881, -0.00963036886324134, 0.0848612242379394, -0.184614356649489, 0.25904172302939, -0.288993376021849, 0.254437908648661, -0.142581730976996, -0.0330984309623288, 0.206364027935802, -0.253198956381471, 0.0769932947759565, 0.188920160250908, -0.154603408595143, -0.161030393082491, 0.0325854306719011, 0.119191267752434, 0.0922856045807742, 0.0571716628074046, -0.130040462992571, 0.198998373790818, -0.256715745325574, 0.28981382087975, -0.277382641964737, 0.195324959134553, -0.0350034362899823, -0.157142688696712, 0.245045000967466, -0.083302979080798, -0.190049206205993, 0.100745510260603, 0.171366647419117, 0.0812027232194978, 0.00991257186300118, -0.0322334954054631, 0.0759636481955999, -0.124257865239343, 0.177880928795453, -0.234459248632815, 0.284734993292882, -0.307505021458093, 0.267308627638937, -0.128895101101712, -0.0879087053719132, 0.227244723802484, -0.071835372107027, -0.201139145536989, 0.0174296066387972, 0.112847449817149, 0.0970196937817123, 0.013358304297849, -0.0325732325128704, 0.0566922581667098, -0.0889591548571636, 0.132956784935785, -0.191152349032791, 0.260672036629152, -0.32414233474661, 0.337069639413799, -0.230987553099094, -0.013500795658052, 0.206566069536364, -0.0297765969855272, -0.206399019462268, -0.0823981869115123, -0.0172696139717228, 0.00246676362832286, -0.00651626938344927, 0.0129080515057394, -0.0238403138170055, 0.0429015797379947, -0.0757267431102765, 0.129967569942016, -0.212023741399845, 0.314617073771157, -0.389196083842969, 0.324811700742985, -0.0449593808712374, -0.189995826664464, -0.054745256625771, 0.150540024381074, 0.100282180697952, -0.00158730786055036, 0.00355704285851566, -0.00522335498295973, 0.00601189929687713, -0.00457496416747483, -0.00211180980942479, 0.0203964304632372, -0.0622133290993153, 0.145240321341511, -0.279953554146627, 0.421258107067016, -0.400409928499652, 0.0646059917093367, 0.16514585654828, 0.164033398578897, 0.00045586808051677, 0.00191124793462794, -0.00454234540335871, 0.00756557213085381, -0.0111861671668544, 0.01553460380335, -0.0203594927612126, 0.024128027806026, -0.0216654387951066, -0.000969072150733292, 0.0733269607046042, -0.232541949947705, 0.44526349399145, -0.452139689026344, 0.0241142725959038, 0.0835707129936931, 0.1684546815832, -0.00114437379058268, 0.00275320151972116, -0.00470117619682267, 0.00725928953888057, -0.0108431207031455, 0.0160781233205372, -0.0237877741950876, 0.0345002861071944, -0.0460732028409267, 0.0457593638429548, 0.00750598545277653, -0.189477338553291, 0.482811004057322, -0.462098108531499, -0.0721595286130255, -0.0759033499259238, 0.000486577014401278, -0.00117167634099952, 0.00200662619616866, -0.00312306100661635, 0.00474989336844293, -0.00731097217703883, 0.0116209525274379, -0.0192241904490457, 0.0326940404289481, -0.0543125672202219, 0.0762163255057659, -0.0458846949782158, -0.171146837040925, 0.556167319176913, -0.369461040840413, -0.111158024129943, -0.000155912497775315, 0.000372998495665334, -0.000631523997537749, 0.000967793304682137, -0.00144730371604558, 0.00220110019611122, -0.00351836687128071, 0.00608721920163425, -0.0116044778629387, 0.0240752697483119, -0.050862887491817, 0.0929017210000717, -0.0833779684963333, -0.209120613369985, 0.65506068482068, -0.0906767085339321, -3.69026491708002e-05, 8.71362219437661e-05, -0.000143888621351305, 0.000211992979372342, -0.000299509657273189, 0.000422214884806931, -0.000618013051270166, 0.000996729872018763, -0.00192162354195972, 0.0046686673428354, -0.0135784652378476, 0.0407247903593205, -0.099620189473613, 0.0981177928809731, 0.362349012413802, -0.589306218727437, 5.7058980433382e-06, -1.32112823495045e-05, 2.09681518516171e-05, -2.88256416233254e-05, 3.61414545464006e-05, -4.10550321098359e-05, 3.90095644000945e-05, -2.05576863605841e-05, -2.54555847368395e-05, 5.34773026054131e-05, 0.000458864688290089, -0.00477927648658511, 0.027196177397784, -0.0950140628513317, 0.0702699546297211, 0.696614504090793).finished();
}
} // namespace mrcpp
} // namespace detail
