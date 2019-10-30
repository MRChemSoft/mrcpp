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
auto get_L_G0_22() noexcept -> Eigen::Matrix<double, 23, 23> {
  return (Eigen::Matrix<double, 23, 23>() << 0.0339717517016027, 0.0588407999692904, 0.0742936264245932, 0.0800037948889865, 0.068906265597178, 0.0318660438635907, -0.0346253730907703, -0.114902269304469, -0.161615833327443, -0.107864276927583, 0.0640047727487355, 0.217950051626652, 0.102133119058666, -0.240190695849874, -0.146237150827974, 0.405311157739791, -0.31530335618963, 0.138384156152883, -0.0392045609971086, 0.00747405747584584, -0.000958375545155701, 7.95788306785348e-05, -3.87659371888857e-06, -2.14557714307727e-54, 0.00462151457776881, 0.0178990489939534, 0.041636525641038, 0.0743293390971467, 0.107857764918552, 0.124342504826107, 0.0981880514693807, 0.0115281199739494, -0.113108577863439, -0.184406895714098, -0.0846367449917118, 0.157465262190681, 0.218886254238054, -0.142694154981673, -0.266312347888965, 0.400808174394942, -0.25732761838054, 0.100105534682779, -0.0259085163211725, 0.00459051433802688, -0.000553303486944985, 4.35403668583563e-05, 0.0313046309711751, 0.0542212113542695, 0.068694948398825, 0.075107721685841, 0.0680535455895309, 0.0399825735291111, -0.0136245225391635, -0.0848280037322188, -0.14264268835971, -0.134423434845161, -0.0213966756080519, 0.146829028385006, 0.192413047718281, -0.0287491387156086, -0.265192777886444, 0.00996953747734386, 0.347307623565237, -0.369402431162008, 0.202524168893423, -0.0707847575947065, 0.0168635954802411, -0.0027908328853353, 0.000317334210234233, -9.10355975202715e-53, -0.00394815360383362, -0.0152911331559164, -0.0356609990334652, -0.0642215204263592, -0.0952413109233903, -0.115508990152277, -0.104725572379329, -0.0443521474699878, 0.0608399196334685, 0.156823527680588, 0.145874990695458, -0.0241723586712417, -0.211236846857001, -0.113052963398928, 0.235223113127254, 0.126189878250867, -0.386981911769356, 0.322516119312817, -0.154574072233629, 0.049062345894549, -0.0108283130435128, 0.00168111899529491, -0.02902600344959, -0.0502745127153594, -0.0638656509679207, -0.0706519397989394, -0.0664443080677417, -0.04499499597134, -0.00159792610007854, 0.0602297571583775, 0.120294753614494, 0.138964931434563, 0.0747403982660197, -0.0673659146255219, -0.183893895431667, -0.113805289448452, 0.140173853781128, 0.218646588899975, -0.145251536236476, -0.2419124646, 0.390767318039772, -0.269583020831342, 0.114902929079825, -0.0334101725634901, 0.00686858402616148, -1.23295326981711e-51, 0.00341208398194494, 0.0132149444379339, 0.0308813630316933, 0.0559963300602338, 0.0844419138927609, 0.106299786503036, 0.105494108723779, 0.0648814883673849, -0.0195400323939716, -0.117944913275781, -0.158625384751437, -0.0690506199313752, 0.119334135808703, 0.200772616388427, -0.0140793800780276, -0.259794353904855, 0.0223305752553413, 0.323991320118172, -0.367851218503088, 0.217344901545401, -0.0834711265473854, 0.0223956292987398, -0.0270619170841706, -0.0468726153399998, -0.0596718416982707, -0.0666269658699973, -0.0644578490520059, -0.0479717647811049, -0.012648249107594, 0.0405372098337363, 0.0985020763679563, 0.131977887589304, 0.103525487018553, -0.00233935575484837, -0.134548749839095, -0.16726668573612, -0.00791705405617799, 0.203251906615983, 0.11865548503284, -0.232188174504104, -0.105991697847448, 0.369007966374286, -0.328174684796753, 0.169969954147839, -0.0594309005051338, 7.20122730734577e-51, 0.00298243369266948, 0.0115509160228768, 0.0270367120006932, 0.0492939947190474, 0.0753176044018885, 0.0975431028523044, 0.103174935035769, 0.0770971494990361, 0.0113000245734692, -0.0794358261516174, -0.145661267219754, -0.118908245526358, 0.0233859185983557, 0.175116919460842, 0.133974878549428, -0.126857505103518, -0.216515740716696, 0.1497991179101, 0.218620136024563, -0.380951837822581, 0.280763784492841, -0.129634503813219, 0.0255216905594211, 0.0442048647439683, 0.0563739605970957, 0.0634166635490391, 0.0627271723948767, 0.0499380464521097, 0.0208414958264745, -0.0250440409431553, -0.0792994396347146, -0.120370926121561, -0.116568190974047, -0.0444344973596506, 0.0773217824762001, 0.163716591157338, 0.10518785948429, -0.0939811921805595, -0.205836777093811, 0.00398161750204866, 0.256826455243301, -0.0347735312090645, -0.305338458887712, 0.369350519619053, -0.233411437844591, 1.94208525563907e-50, -0.00268581350212819, -0.0104021109647615, -0.0243805838064101, -0.0446519904600921, -0.0689572825014917, -0.0913388553962737, -0.101347064993389, -0.0856849541444948, -0.0346069615591213, 0.0462248277847075, 0.12367700856024, 0.140073100121304, 0.050792764489282, -0.10690510357101, -0.182410631002471, -0.0360143157996741, 0.199875046942233, 0.127511982568985, -0.235554047893049, -0.0964174200123429, 0.370666422044367, -0.348382325289261, -0.0253305758808503, -0.0438738444106115, -0.0560360142687324, -0.0634399137536062, -0.0639211726115269, -0.0536290169604623, -0.0285468300548957, 0.0127298163464909, 0.0647825550115645, 0.1112740703769, 0.125253759653548, 0.0802998942969501, -0.0243068478733527, -0.136156236410413, -0.156361592894118, -0.0214918405690899, 0.16961876282998, 0.164341803455418, -0.114497576695644, -0.238346309794226, 0.149437173623672, 0.249364291238798, -0.435158518821773, -2.57872214732388e-50, -0.00258677060121636, -0.0100185194589699, -0.0235106557166035, -0.0432366650759452, -0.0674202079860697, -0.0911037527089616, -0.105278792392789, -0.0977428810307564, -0.0575639074613643, 0.0157241018549903, 0.100734120527409, 0.150337101206697, 0.111221138724752, -0.0249511233763492, -0.168228531018999, -0.158542746010715, 0.0572401277143087, 0.237148516257793, 0.0457676643454246, -0.291257701198554, -0.0366947544745407, 0.462018638097722, 0.0268571669899428, 0.0465179777739426, 0.0595015484812002, 0.0677863689322626, 0.0695249156786556, 0.061171717194508, 0.0387442188023426, 5.18915270971933e-05, -0.0517935145223517, -0.104474271623551, -0.135056317149748, -0.116019603101349, -0.0330109170721251, 0.0893781097875767, 0.17470669730686, 0.130395395395521, -0.0534436223201047, -0.21664937678209, -0.126293526601918, 0.183774138351743, 0.24975483138048, -0.185835819263129, -0.424937865546294, 1.9228333814119e-50, -0.00256110284035796, -0.00991910864863089, -0.0233074923050542, -0.0430467571980886, -0.0677936377274388, -0.0934690504744307, -0.112351605992408, -0.113325069807829, -0.0847458001348089, -0.0211822645745976, 0.0667273440144551, 0.144684855317087, 0.160586679363843, 0.0753588088686518, -0.0861254035479707, -0.206452110989289, -0.140907367514534, 0.107591382597995, 0.272710209344251, 0.0613696647348518, -0.340259712448485, -0.344720844347002, -0.0290274444257063, -0.0502770085592053, -0.0644114458589876, -0.0738656560619041, -0.0771622050662393, -0.0711219440145972, -0.0518888271146119, -0.0166243798703916, 0.0337346200753068, 0.0910405236546537, 0.13751816498465, 0.147710193509878, 0.0987642121325688, -0.0100330420894497, -0.136966839047803, -0.198790766615154, -0.11827116779154, 0.0869459556585869, 0.251073524360979, 0.164946801733916, -0.162376646829999, -0.392669188558988, -0.249813013671692, -9.04834034097589e-51, -0.00250174404499391, -0.00968921302273499, -0.0227986299297911, -0.0422977856908693, -0.0673094604407305, -0.0947365146644201, -0.118400724348983, -0.128835075324337, -0.114745493546849, -0.0672408955603447, 0.012898930879666, 0.107065976352818, 0.175835690365908, 0.170262917503264, 0.0631766232098569, -0.109943516126159, -0.236891118747066, -0.187883451938087, 0.0553167052438364, 0.316737545203193, 0.356231482387661, 0.16267526762291, -0.0317020307884055, -0.0549095280286311, -0.070464073285962, -0.0813683213145946, -0.0866163749425192, -0.0835351672003269, -0.06859516181059, -0.0385791819448457, 0.00747728705535096, 0.06562561616538, 0.124397443792952, 0.164072142149997, 0.160542292389621, 0.0960769858994165, -0.0241401562214765, -0.157950192749876, -0.230516679051832, -0.171439140127003, 0.0231730457863305, 0.251740260278744, 0.357745007831904, 0.270088490969376, 0.0951754518330247, -2.75277262465619e-51, 0.00239163031747944, 0.00926274438988264, 0.0218267984135899, 0.0406872344278068, 0.0654476102863204, 0.0940698396363644, 0.122152956294595, 0.142479217735925, 0.14541240857738, 0.120896605682111, 0.0626252568341059, -0.0259345105624669, -0.125794978581951, -0.201425325429203, -0.209978319988907, -0.12311827296245, 0.0458824000157681, 0.230072526221021, 0.335278814175329, 0.304852624444435, 0.174811398504645, 0.0497695706493358, -0.0354062352178421, -0.061325398302037, -0.0788359858284245, -0.0916956089504023, -0.09950362385378, -0.100268346233545, -0.0910540353641963, -0.0686806982485718, -0.030861398666445, 0.0221174395389338, 0.0853758608548884, 0.147763430717811, 0.191780352173573, 0.196677126529909, 0.145771263633485, 0.0369514172173007, -0.108217788039637, -0.244608126215466, -0.319584641466451, -0.302410823109092, -0.209053344548161, -0.0969267608024939, -0.023019479898837, -5.06854621312539e-52, -0.00220470526511124, -0.00853878677507163, -0.0201515386450579, -0.0377507102539996, -0.0614025446155277, -0.0901499675090595, -0.121526781578211, -0.151105598784215, -0.172320994903869, -0.17691484538644, -0.156380685046753, -0.104616748239766, -0.0215053431461594, 0.0837030195320807, 0.191243620843068, 0.274274042605144, 0.307924299101683, 0.28153234318765, 0.207289261991828, 0.117120871959358, 0.0456122990603728, 0.00925237624279761, 0.0416291965020901, 0.0721038834198886, 0.0928632094041763, 0.108824150853912, 0.120418010936839, 0.126599330202789, 0.125425705616724, 0.114409111224535, 0.0909810808934999, 0.0531899670966012, 0.000642703377144758, -0.0644327973386674, -0.13631440099798, -0.205601850136131, -0.260298449796814, -0.288517410689149, -0.282511391450768, -0.242687653726193, -0.179361351550131, -0.110151431295044, -0.052945562479201, -0.0177736110467187, -0.00313603065909749, 4.6760740640664e-53, -0.00185268918054704, -0.00717543434195734, -0.0169610621725911, -0.031937800117643, -0.0525445479690375, -0.0788144215968882, -0.110199457569754, -0.145369799959397, -0.182041351812567, -0.216913994140171, -0.245825410906204, -0.264223328384918, -0.268005483047355, -0.254650442930867, -0.224369377000182, -0.180811360445574, -0.130786220840072, -0.0826949351109797, -0.0439511048106904, -0.0184387225718671, -0.00543695848582089, -0.00084662663921801, -0.0598848615020997, -0.103723622725862, -0.133844338992487, -0.158071984667376, -0.178402399777739, -0.195394471124694, -0.208953907165694, -0.218599173644916, -0.223614681704123, -0.223196715309449, -0.216629121623011, -0.203495737707548, -0.183910376814418, -0.158716788894623, -0.12958551595344, -0.0989257308542845, -0.0695565044670532, -0.0441565855825184, -0.0246252486463892, -0.011591279992179, -0.00432145333548557, -0.00113549948407031, -0.000157902675727757).finished();
}
} // namespace mrcpp
} // namespace detail
