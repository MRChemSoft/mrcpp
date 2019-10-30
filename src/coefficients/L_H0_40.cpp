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
auto get_L_H0_40() noexcept -> Eigen::Matrix<double, 41, 41> {
  return (Eigen::Matrix<double, 41, 41>() << 0.707106781186548, 1.70550287454142e-309, -3.38949490201654e-308, -6.14595630465374e-310, -3.06683219602222e-308, 3.0729781523269e-311, -2.49833123784175e-308, -3.07297815232685e-310, 2.61203142947784e-308, 1.22919126093076e-310, -1.27036916817193e-307, -1.84378689139614e-310, 9.35721847383533e-307, 1.22919126093076e-310, 1.33009557154055e-305, 6.14595630465379e-311, -5.158848116607e-305, -6.14595630465379e-311, 1.88514118758325e-304, -1.22919126093076e-310, -3.50428968847049e-304, -6.14595630465379e-311, -1.44755585924384e-304, 6.14595630465379e-311, 7.23230896489109e-303, -1.22919126093076e-310, -9.08886716270722e-302, -1.22919126093076e-310, -1.3242379380871e-300, -6.14595630465379e-311, 3.24123220857491e-300, -1.22919126093076e-310, -7.67529519386845e-299, 6.14595630465379e-311, -4.27849236370133e-299, -4.6094672284901e-311, -2.60200385388945e-298, 3.0729781523269e-311, 4.68543214496436e-297, 0, 1.46029301998771e-296, -0.612372435695794, 0.353553390593274, 2.93394389666544e-308, -2.58619920688484e-308, 2.60238635937864e-308, -1.74382507278869e-308, 2.09096356961346e-308, 1.9415099974043e-308, -2.27939955021982e-308, 3.91582643734887e-309, 1.09439015271236e-307, -1.92764648454267e-306, -8.1012388859763e-307, -4.03548172453719e-307, -1.15187808956006e-305, -1.07551566882053e-305, 4.4676595500273e-305, 4.46378417646641e-305, -1.63258318204848e-304, 6.31219748568653e-304, 3.0347949583257e-304, -2.49278543873169e-303, 1.25362858191595e-304, 1.00841019557782e-303, -6.26336283474557e-303, -4.93143288949391e-302, 7.871189900902e-302, -2.06494841539194e-304, 1.14682369613251e-300, -1.39850291343698e-300, -2.80698943266615e-300, -1.7097289264258e-300, 6.64700061938055e-299, -1.21307441529213e-298, 3.70528307691584e-299, -7.41424854367929e-299, 2.25340143820119e-298, 4.29034799916118e-298, -4.05770326524743e-297, 1.83418473749138e-296, -1.26465085227851e-296, -1.16312223065572e-308, -0.684653196881458, 0.176776695296637, 4.79077293947759e-308, -1.97285197379385e-308, 2.9162562665582e-308, -2.82713990014072e-309, -3.93955799128305e-308, 1.46888355681224e-308, -7.62098581777065e-309, -5.76490701376521e-308, 3.73409867201848e-306, 3.74718955894739e-307, 7.82687535397654e-307, 4.99027068112665e-306, 2.08266792104541e-305, -2.15671440260387e-305, -8.64426910462658e-305, 7.5811723128291e-305, -1.2223497667511e-303, -5.86541183721522e-305, 4.82725767156541e-303, -2.18103272251013e-304, -1.9527790346646e-303, 3.93343034992818e-303, 9.54967880665561e-302, -4.68248220347239e-302, 3.99876137377015e-304, -6.25862746284626e-301, 2.70818924688883e-300, 1.715403466354e-300, 3.31087583031952e-300, -3.98909046777308e-299, 2.34910850406474e-298, -1.49770325581551e-299, 1.43576305671146e-298, -1.33971196344236e-298, -8.30822317509602e-298, 2.3475412662683e-297, -3.55188347108609e-296, 7.99379708585455e-297, 0.233853586673371, 0.405046293650491, -0.522912516583797, 0.0883883476483184, 3.88545677045283e-308, -2.25494176512073e-308, 4.36602974048566e-309, 2.46740939518396e-308, -3.20876537443282e-308, 1.62027574609797e-308, 1.24582135441604e-307, -2.65239547766149e-306, -8.00759047697518e-307, -3.47909862356287e-307, -1.03633489728284e-305, -1.46052203127505e-305, 4.67365662846471e-305, 6.98852967816932e-305, -1.61910448287598e-304, 9.15808776581586e-304, 5.76060330786288e-305, -3.69767145743892e-303, 5.97284863916672e-304, 1.29711499293145e-303, -9.24338420840252e-303, -7.18202957777162e-302, 1.08451086476184e-301, -2.4706350244503e-302, 1.4133764374319e-300, -1.85677772202844e-300, -4.00229508763904e-300, -2.99590651768728e-300, 9.26152344432761e-299, -1.75767832609758e-298, 3.01528762849618e-299, -9.81894789158906e-299, 3.1023895425952e-298, 6.39002484648587e-298, -5.39455968262358e-297, 2.66181733844103e-296, -1.88165054496976e-296, -2.93162115731984e-308, 0.153093108923949, 0.592927061281571, -0.350780380010057, 0.0441941738241592, 1.13700191636094e-308, -1.45659164420294e-308, 1.29065082397729e-309, 5.67271766919541e-308, -4.20383411238316e-308, -2.01833205044829e-307, 9.23245556085086e-307, 1.3108710202196e-306, -6.2799381520952e-307, 1.73144495610336e-305, 4.41666857921332e-306, -7.61702782191045e-305, -5.50699195744042e-305, 2.6655854489297e-304, -4.91251204196018e-304, -1.5388687904446e-304, 2.26149091152386e-303, -8.61638735925489e-304, -1.26241691530934e-304, 1.42989424304902e-302, 3.94600583081543e-302, -1.70303654406886e-301, 9.7022009515188e-302, -2.2400754418326e-300, 4.0480860054525e-301, 6.26269925056504e-300, 3.37579045806232e-300, -1.45079676257083e-298, 9.34888356718162e-299, -5.00957095203364e-299, 2.04740775070318e-299, -4.85452571698743e-298, -3.99521036240262e-298, 8.4730590424796e-297, -1.43016924655913e-296, 2.93716124805139e-296, -0.146575492494482, -0.253876200144874, -0.163876382526586, 0.581703452155821, -0.219863238741723, 0.0220970869120796, 8.21205395925726e-309, 5.63387330473865e-309, -6.03391264030289e-308, 7.00578999625989e-308, 1.45240471147039e-307, -1.59420008601697e-306, -7.90167596357971e-307, 1.11450651590382e-306, -1.027502557882e-305, -7.36589958188092e-306, 4.97340660916424e-305, 9.72619811486247e-305, -1.70430898651777e-304, 8.58655759301333e-304, -8.61325271987532e-305, -3.94996134274898e-303, 8.79424177967012e-304, 1.84779397088001e-304, -1.06469217529219e-302, -6.88760746710166e-302, 1.28158113648515e-301, -1.76206691114704e-301, 1.55549962551298e-300, -6.45476734848901e-301, -4.79464056693276e-300, -6.0389540040921e-300, 1.09014730537043e-298, -1.62494628358778e-298, 2.20624811267646e-299, -3.16833371843158e-299, 3.57724151045169e-298, 7.06225791195403e-298, -6.12828230802013e-297, 2.48365953589263e-296, -2.31713347076265e-296, -6.06144940546476e-308, -0.0689981317681863, -0.26722861525761, -0.421585548851001, 0.478033079399324, -0.132121363478811, 0.0110485434560398, -4.08706094259474e-308, 7.19691483274953e-308, -6.31804308118405e-308, -5.61125810614887e-308, 2.80501445744397e-306, -2.6304692983918e-308, -8.82743704037417e-307, -1.24160609266615e-306, 1.33103590285517e-305, -9.41462170572079e-306, -1.33637182211887e-304, 1.97797770135193e-305, -1.30521889150402e-303, 5.72562390485276e-304, 5.69576850853295e-303, -1.11925933948373e-303, -7.26420936671171e-304, 6.21822519127022e-303, 1.03275138821504e-301, -7.57187931959865e-302, 2.08317495076322e-301, -5.99713899080406e-301, 1.3328544851183e-300, 3.03895817765409e-300, 7.91646333987729e-300, -6.43076569569816e-299, 2.43773781993559e-298, 2.54200710372217e-299, 6.40939951916852e-299, -1.94929872944953e-298, -1.04842883935927e-297, 3.04225041154507e-297, -3.70095766688913e-296, 1.63202870881394e-296, 0.106977062012728, 0.185289706650491, 0.181798066847189, -0.0566069404148025, -0.534885310063639, 0.3548027758708, -0.0771422564770762, 0.0055242717280199, -8.77481709104394e-308, 4.70026412983484e-308, 2.85542090219886e-308, -1.50444847739707e-306, -1.12473881292181e-307, 9.07519590390553e-307, -1.39275052308835e-307, -4.8772772744656e-306, 2.05031157377388e-305, 1.08356686917886e-304, -5.32737182298497e-305, 9.23458210974755e-304, -7.74794127185763e-304, -4.00696596807494e-303, 1.61988686791509e-303, 3.76597294584235e-305, -9.66238447121417e-303, -7.19583508973112e-302, 1.20776641965196e-301, -2.35434705391649e-301, 1.01521481452692e-300, -1.27106441575613e-301, -4.78471472093222e-300, -7.44441558351587e-300, 1.01929052136457e-298, -1.60886123184507e-298, -3.24817048439462e-299, 8.41338495395248e-300, 3.08470355671726e-298, 8.46284938545003e-298, -4.89665417931219e-297, 2.41419350385842e-296, -2.55618760122076e-296, -4.45735480995013e-308, 0.0394511911654769, 0.152793806371937, 0.301313449619951, 0.204994402553155, -0.528802957971469, 0.246372609862836, -0.0441077726196727, 0.00276213586400995, -1.84378689139612e-308, -2.02201962423108e-308, -5.34636738941829e-307, 5.21914609391196e-307, -1.45161341959617e-306, 7.83197649770941e-306, -8.47674877363063e-306, -4.98165405038753e-305, -8.59921309895716e-305, 1.58419030142638e-304, -4.01315090126624e-304, 6.09146748524795e-304, 1.89601627288564e-303, -1.73213072190783e-303, 1.22802709388752e-303, 1.30796385996837e-302, 3.02224743155246e-302, -1.73990374823645e-301, 3.23831546422723e-301, -1.68948225221895e-300, -1.86667180948831e-300, 6.67405802311987e-300, 7.89530497969186e-300, -1.44869862902802e-298, 4.81928683563569e-299, 1.66166455434959e-299, -1.26513122411107e-298, -4.39516835510032e-298, -6.09181207998631e-298, 7.27395366551009e-297, -6.74003262767907e-297, 3.49794666389029e-296, -0.0842790976968415, -0.145975679226991, -0.161531821313557, -0.0637090094933597, 0.216717679791878, 0.39931734961447, -0.455808912293701, 0.163205770906627, -0.0248208301311305, 0.00138106793200498, -8.01447106711937e-308, 1.00443651978658e-307, 4.3873005111338e-308, 1.34334519699715e-306, -4.75518737232081e-306, 7.36696768186527e-306, 3.30526135383726e-305, 1.22973898087417e-304, -9.96685811077025e-305, 6.76767409975831e-304, -5.80672696217298e-304, -2.95513574171096e-303, 1.38271401991302e-303, -1.11930354523475e-303, -7.60925645174419e-303, -5.18719981363188e-302, 1.27207828384125e-301, -4.11001761226456e-301, 9.72359865754583e-301, 2.06788632346691e-300, -4.87885666585693e-300, -1.05092158958184e-299, 1.01248646496331e-298, -9.01593052938879e-299, -3.69506269185855e-299, 1.51572618642185e-298, 2.70859407975614e-298, 9.33657969059372e-298, -4.41946560045054e-297, 1.26654317281494e-296, -2.77972641456424e-296, 7.69166431527416e-308, -0.0255777360424, -0.0990621457259044, -0.210503024770429, -0.253153933531557, 0.0104952076465843, 0.490608224764219, -0.359502355151361, 0.10437804074912, -0.0137934051677984, 0.000690533966002488, 1.31437421531325e-306, -1.05065123028056e-306, 9.79665434961807e-308, -4.05374985942352e-306, 1.2407456587835e-306, 2.11390167098566e-306, -1.5301820958036e-304, -3.37579556541347e-305, -1.11711735079243e-303, 8.70737086719772e-304, 4.27389912052835e-303, -1.33369605712251e-303, -1.21373110784639e-304, 1.31105650518087e-303, 8.6295591085309e-302, -6.77754473605695e-302, 3.93854793409416e-301, 8.41742018572169e-302, -1.19847324909818e-300, 2.69255465570798e-300, 1.14160121805171e-299, -4.65427777968871e-299, 1.63898449670884e-298, 7.98338911205761e-299, -1.24387795015792e-298, -5.59931616375733e-299, -1.34289776149968e-297, 6.25888574971651e-298, -2.30092084632554e-296, 1.9868298777535e-296, 0.0695453758035306, 0.120456124323186, 0.1406978842158, 0.0963808314913037, -0.0645778489604212, -0.291065800584817, -0.208957848692535, 0.493804817536819, -0.266261292339, 0.0649588425253176, -0.00758802259175784, 0.000345266983001244, 2.06808548734581e-307, -6.15335545986114e-307, 9.65962353166611e-307, 1.97683340111246e-306, 8.50681834500399e-306, 1.19209508013415e-304, 1.05034119559416e-305, 7.81192874835825e-304, -9.8290365483183e-304, -2.21626086245511e-303, 1.69562707371887e-303, 6.04681247441824e-304, -3.7566989210497e-303, -5.96509468020658e-302, 1.06158132797061e-301, -3.30459518855338e-301, 2.32317104720435e-301, 1.82957777054404e-300, -3.98516826084512e-300, -9.29554276877786e-300, 7.48473777183738e-299, -8.37060894152159e-299, -8.01748554011347e-299, 1.84919835049871e-298, 1.23472629781877e-298, 1.25124577749629e-297, -1.87540922660895e-297, 9.63210053604071e-297, -2.80960880279682e-296, 2.60182914201212e-306, 0.0179405987020252, 0.0694836399939341, 0.153009670967929, 0.220107701663293, 0.143140129389143, -0.180507282585682, -0.371745290494229, 0.441307394385264, -0.18812298593045, 0.0395553021143058, -0.00413960570259111, 0.000172633491500622, -1.41080426973327e-306, 1.18681489221016e-305, -1.30301034210595e-305, -4.13410823810686e-305, -8.93078744159323e-305, 9.55182403426588e-305, -2.52827249317113e-304, 6.3731059329087e-304, -3.19038312642341e-304, -1.63600685903091e-303, -3.25558742064637e-304, 7.38128344252081e-303, 2.00933850664687e-302, -1.59012514290277e-301, 3.43865280619626e-301, -9.85390233741296e-301, -3.4307138262331e-300, 5.5252451938987e-300, 8.35051115739113e-300, -1.13643112028911e-298, -3.36007061249571e-299, 4.06186782318304e-299, -3.00939560521469e-298, -2.32580939337756e-298, -1.11849601574526e-297, 4.1815964560812e-297, 9.37664113407602e-297, 3.67703379294345e-296, -0.0592039757166833, -0.102544293951369, -0.123357924633427, -0.103239338518338, -0.00642871453489237, 0.168093812205297, 0.274914431263867, 0.00810642466595898, -0.454921548126905, 0.363635843788664, -0.128117052855236, 0.0236610244206259, -0.00224257483775315, 8.6316745750311e-05, -1.41348669156855e-305, 1.39886767022282e-305, 2.54190683628628e-305, 1.31443175682915e-304, -6.36999235303254e-305, 4.78937385573859e-304, -3.43137739513277e-304, -4.67423371935713e-304, 1.06363343509513e-303, 1.14829816292217e-303, -2.54239560306446e-303, -4.29637463075026e-302, 1.1551664593916e-301, -3.71931191156972e-301, 4.48598578352686e-301, 3.3049575039304e-300, -3.39922929728526e-300, -9.86481436979621e-300, 6.7997365848117e-299, -3.22950566442564e-300, -2.57543757811107e-299, 3.22065946986574e-298, 5.2090120170763e-299, 1.52296520498296e-297, -1.51871432842022e-297, -4.90186887249095e-297, -2.89195197583771e-296, 2.12689729559185e-305, -0.013284298602997, -0.0514498672554539, -0.115603486198934, -0.181283661577571, -0.176520799885996, -0.00586588815455456, 0.273642327022754, 0.199320449302147, -0.467026783576276, 0.282240662668031, -0.0847113439910263, 0.0139448970444285, -0.00120766350936337, 4.31583728751555e-05, 6.02746226710002e-306, 2.21174529535574e-306, -1.52770650460409e-304, -7.71158335965754e-305, -8.98126971208403e-304, 5.26216594425764e-304, 1.83622146502335e-303, -9.91220140112371e-304, -3.29026548741888e-303, -3.08410473016125e-303, 8.37577009694954e-302, -5.63982487981028e-302, 2.65991165448896e-301, 5.34268367191578e-301, -1.67138257450711e-300, 7.76329618826995e-301, 9.2460687752745e-300, -8.57066949982246e-300, 8.77385723247737e-299, 4.02589468538712e-299, -2.62216397846762e-298, 1.80348677540867e-298, -2.0202184420033e-297, -2.27486831797866e-297, -6.51173474398915e-297, 2.01511919970059e-296, 0.0515434170219584, 0.0892758170777425, 0.109344092753528, 0.101404123890135, 0.0407301127516175, -0.0856743299651528, -0.220498364544324, -0.184736218894801, 0.159017677884524, 0.344078934089262, -0.430008059376046, 0.209163768146521, -0.0546672604778347, 0.00811598337070386, -0.000647015840060858, 2.15791864375777e-05, -1.26119920930425e-306, 1.10291572246608e-304, 9.52045651374544e-305, 5.69982249846879e-304, -7.57200883602471e-304, 7.39835143155158e-306, 1.66819934262338e-303, 4.79079382324506e-303, -8.44403076563975e-304, -6.63499704873149e-302, 9.36315024385086e-302, -8.9995322038301e-302, -2.22838651641182e-301, 1.26394656950687e-300, -1.61982451230618e-300, -5.25948153439036e-300, 3.20603207222426e-299, -2.69885999902547e-299, -1.91053692199637e-299, 2.72695629458144e-298, -1.19912028176848e-298, 2.03305328323856e-297, 9.91882718242649e-298, -4.56265767580689e-297, -2.75334532214823e-296, -4.75843138091343e-305, 0.0102345192513237, 0.0396381226168159, 0.090193190712834, 0.148631623034035, 0.172042866302963, 0.0907318468848569, -0.119866611828528, -0.278556309580142, -0.0124278691962095, 0.424348123789452, -0.366809959680973, 0.149375400964147, -0.034568341247478, 0.00467292785582249, -0.000345098354428059, 1.07895932187889e-05, -6.88063162939944e-305, -1.24092389151633e-305, -2.61232643538046e-305, 5.8482922759286e-304, -2.40195201048037e-303, -1.96900472806626e-303, -5.00276946662518e-303, 7.02964925164242e-303, 3.13087929886236e-302, -1.49855939164346e-301, 2.529390012406e-302, -6.2888510882994e-301, -2.17106535595784e-300, 2.84913446832408e-300, 3.06485824293889e-300, -7.06619132898117e-299, -8.07684799689046e-299, -4.93646052849407e-299, -3.508914976931e-298, -1.01893798756588e-299, -1.9739679779191e-297, 1.80798406343793e-297, 2.26847658115431e-296, 3.55873389642499e-296, -0.0456399474313251, -0.0790507078058277, -0.0979718639481214, -0.0966015606026778, -0.057448884878591, 0.0338731370251748, 0.156041624778151, 0.210533095503393, 0.0531636663628142, -0.258761013058376, -0.183866825334376, 0.444497812568922, -0.295414932464577, 0.103484562945811, -0.021484146740091, 0.0026655154492612, -0.000183343732310454, 5.39479660939444e-06, 2.06797207524587e-305, 2.19865276630413e-304, -3.97729574749798e-304, 1.39306134138418e-303, 1.58055377845344e-303, 5.60022175095304e-303, -4.61299078110391e-303, -5.42052757366519e-302, 1.12360008120684e-301, -2.18081416870279e-302, 2.65154867455442e-301, 1.69586039603402e-300, -5.08952476428691e-301, -4.17084290517521e-300, 2.69338754674981e-299, 4.02417929864994e-299, 9.78213376242113e-299, 3.40020767414457e-298, -1.38232175202296e-298, 2.39261541981614e-297, 1.60435128762113e-298, -1.7861633321012e-296, -2.71694291627202e-296, 4.46130742809856e-304, -0.00812779240754105, -0.0314788046358376, -0.072235100975407, -0.122860358428221, -0.156506062959549, -0.124232982942462, 0.0142344139090215, 0.201327606218267, 0.205689264120524, -0.147671631501979, -0.317271004569918, 0.420035413429816, -0.227379255429147, 0.0698910274057256, -0.0131549524088687, 0.00150807484138255, -9.70688678741415e-05, 2.69739830469722e-06, -3.7247568184354e-304, 4.81558232129259e-304, 6.16205194962e-304, -1.68663970545117e-303, -7.45979828015101e-303, 5.36816340885712e-304, 9.67999432688522e-302, -5.49251192459179e-302, -1.16870092564158e-301, 6.62407076601897e-301, 4.64806171311401e-301, -2.32337860151783e-300, 3.18154055954437e-300, 3.33348800703785e-299, 5.90768003587565e-299, -1.02715449737296e-298, -2.27730320729432e-298, 3.481713085227e-298, -2.88084738190726e-297, -3.56613797435336e-297, 4.20696464336786e-297, 1.76373415622515e-296, 0.0409507954903469, 0.0709288583996433, 0.0886307272903126, 0.0909639704527409, 0.065258486931141, -0.00145260345997924, -0.103223850969413, -0.186198844713189, -0.144326934607074, 0.080035135058834, 0.275687717919289, 0.00888656633110127, -0.396916132760715, 0.368867993548947, -0.168683485161818, 0.0461917552579768, -0.00795115139099776, 0.000847091515042966, -5.1232818674329e-05, 1.34869915234861e-06, -2.77900567305001e-304, -1.27255886676965e-304, 3.42285807120714e-303, 8.17863767652022e-303, -5.49188887577461e-303, -8.24685684550476e-302, 8.92225696717432e-302, 3.38843985298e-301, -4.08087022256305e-301, -1.54815876550106e-300, 1.59450283417281e-300, 1.40016749678047e-300, -1.1226331637726e-299, -2.01003590412706e-299, 1.26063021133781e-298, 1.73998536507039e-298, -2.61705956117014e-298, 2.89200202421548e-297, 2.03944943222768e-297, -1.21066825266345e-296, -2.41988911227385e-296, 1.82945959170762e-304, 0.00661134635589627, 0.0256056343323954, 0.0591088434232506, 0.102723890512604, 0.138970988827876, 0.133391498453093, 0.0478956046615338, -0.110238279264893, -0.219985482132738, -0.0834232810706853, 0.247772805434256, 0.165701105440151, -0.424277917160628, 0.306332923264078, -0.121360459018218, 0.0299640781761113, -0.00475142615106978, 0.000472770516148696, -2.69655523597692e-05, 6.74349576174305e-07, -3.06503007871457e-303, -3.05776254684388e-303, -8.52104483180298e-303, 1.41322151122363e-302, 4.43115346122147e-302, -1.43263632721607e-301, -4.11382352530284e-301, -4.63981743298652e-301, 8.87010590468153e-301, -2.58104468091201e-301, -3.63377760249042e-300, -3.03611605565486e-299, -7.69783124039556e-299, -1.9504924212649e-298, -2.07272085251622e-298, 8.55580395585388e-299, -2.80387702185025e-297, 1.30579421149053e-297, 2.87789900084565e-296, 3.15911981097937e-296, -0.0371360054612984, -0.0643214482491242, -0.0808534054499439, -0.085324656829441, -0.068359789048132, -0.0189633117232361, 0.0635837506744712, 0.149958365846225, 0.168394422018244, 0.0429448711120381, -0.178415943868476, -0.215433292864517, 0.142996615456625, 0.291388108825173, -0.410468654023286, 0.243326247690153, -0.0850729385293069, 0.0191236468801429, -0.00281084222623789, 0.000262346470332942, -1.41573265450511e-05, 3.37174788087152e-07, 1.5102807047194e-303, 8.72717959838757e-303, -1.28846373068885e-302, -6.02998206280216e-302, 1.06892953949985e-301, 3.61025260393191e-301, 2.18794253601394e-301, -1.11475786808329e-300, 2.32783598844808e-300, 1.70584765204416e-300, -7.96838230615545e-300, 3.49062170116532e-299, 2.44728549786843e-298, 1.68985287522679e-298, -1.86765969160996e-298, 3.13122508186261e-297, 1.9777131680482e-299, -2.35879686988478e-296, -2.24989423497758e-296, -7.1567566933082e-304, -0.00548335795155571, -0.0212369540276693, -0.0492384817480257, -0.0868997867372054, -0.122464688632806, -0.131322704818521, -0.0816110359705669, 0.0384968769960545, 0.172066290983573, 0.175488770814139, -0.0501366381202787, -0.270683109786265, -0.000803764424668383, 0.371470721264963, -0.369741542142612, 0.186340754749332, -0.0583166494034511, 0.0120315623487937, -0.00164793638714457, 0.000144828248486399, -7.41592932462492e-06, 1.68587394043576e-07, -1.36196950887332e-302, 1.15074676377045e-302, 9.48782432904076e-302, -4.35553633750766e-302, -4.35812407765944e-301, 6.64790191287033e-301, 3.08677385773991e-300, -4.84098762311282e-300, -1.86024617382103e-300, 6.53147855615526e-299, 7.40009770602419e-299, -2.3511571391679e-298, -2.08967942457617e-299, 3.71748896877852e-298, -3.5016771228283e-297, -3.19582640677667e-297, 8.09104528168185e-297, 1.2025053268314e-296, 0.0339717517016025, 0.05884079996929, 0.0742936264245927, 0.0800037948889859, 0.0689062655971775, 0.0318660438635904, -0.0346253730907701, -0.114902269304468, -0.161615833327441, -0.107864276927582, 0.0640047727487351, 0.21795005162665, 0.102133119058666, -0.240190695849873, -0.146237150827973, 0.405311157739788, -0.315303356189628, 0.138384156152882, -0.0392045609971083, 0.00747405747584578, -0.000958375545155694, 7.95788306785343e-05, -3.87659371888854e-06, 8.42936970217881e-08, -1.70897675169494e-302, -6.87456493787737e-302, 5.02837975182181e-302, 6.35449051703539e-301, -4.90961459649724e-301, -4.12987979063602e-300, 3.67692820886795e-300, 6.16910797786732e-300, -4.11088710047928e-299, -5.13489403860863e-299, 2.24388852701092e-298, -7.20330639558236e-299, -2.69288662172211e-298, 3.41604425061132e-297, 1.71333270812724e-297, -1.3228653180322e-296, -1.75047740348496e-296, 1.32552653855859e-302, 0.00462151457774989, 0.0178990489938801, 0.0416365256408676, 0.0743293390968425, 0.10785776491811, 0.124342504825598, 0.0981880514689788, 0.0115281199739022, -0.113108577862976, -0.184406895713343, -0.0846367449913654, 0.157465262190036, 0.218886254237158, -0.142694154981089, -0.266312347887875, 0.400808174393301, -0.257327618379487, 0.100105534682369, -0.0259085163210665, 0.00459051433800809, -0.00055330348694272, 4.35403668581781e-05, -2.02260965120282e-06, 4.2146848510894e-08, -1.30184327417134e-302, -5.54864418500537e-302, -6.69963370913309e-301, -3.88000206905987e-301, 3.11887838191919e-300, -1.79268318848234e-300, -7.98372820964385e-300, -6.55564219160021e-300, -3.98684596132857e-299, -2.57529760351047e-298, 5.13557826296674e-299, 6.5777429192356e-299, -3.22863515777757e-297, 1.81167539397291e-297, 2.87976934364531e-296, 2.38938861079154e-296, -0.0313046309536483, -0.0542212113239122, -0.0686949483603631, -0.0751077216437837, -0.0680535455514093, -0.0399825735066795, 0.0136245225316148, 0.0848280036848168, 0.142642688279874, 0.134423434769748, 0.0213966755957234, -0.146829028303034, -0.192413047610126, 0.0287491387003486, 0.265192777737401, -0.00996953747334031, -0.347307623368017, 0.36940243095303, -0.202524168778998, 0.0707847575547383, -0.0168635954707225, 0.00279083288376033, -0.000317334210055171, 2.37299576748381e-05, -1.05346045745216e-06, 2.1073424255447e-08, 1.61276888115294e-301, 3.29471992296681e-301, 4.45269073779167e-301, -2.64363529937217e-300, 3.27929976511564e-300, 5.13901116198157e-300, -2.36391935242816e-299, 2.83996096303603e-300, 2.77675969404067e-298, -7.52766773541279e-299, -1.44438042947437e-298, 3.41025858894103e-297, -7.74076887345373e-298, -2.37525442819545e-296, -1.37326317676783e-296, -9.88951117577234e-302, -0.00394815347735416, -0.0152911326660635, -0.0356609978909838, -0.0642215183683974, -0.095241307869668, -0.115508986443972, -0.104725569006555, -0.0443521460195627, 0.0608399177151585, 0.156823522657345, 0.145874985959225, -0.0241723579915451, -0.211236840069211, -0.113052959575701, 0.235223105585243, 0.12618987380636, -0.386981898851779, 0.322516108630937, -0.154574067130241, 0.0490623442773783, -0.0108283126869318, 0.001681118939966, -0.00018090493786476, 1.2887014923139e-05, -5.47807706733102e-07, 1.05367121277235e-08, -2.71585778154287e-301, 3.02194434248811e-301, 3.9623096640592e-300, -5.37147540457151e-300, -4.28182787976618e-300, 7.533123888746e-299, 1.07408488270995e-298, -2.32434461040218e-298, 2.10401023231998e-298, 3.34025630641353e-298, -3.61486782908074e-297, -2.49727176389039e-297, 7.06080897178395e-297, 1.90648673628981e-297, 0.0290259753654005, 0.0502744640721163, 0.0638655891657856, 0.070651871387702, 0.0664442436081886, 0.0449949520307467, 0.00159792382538831, -0.0602296998266835, -0.120294637837257, -0.138964796248224, -0.0747403231416156, 0.067365852973061, 0.183893717240485, 0.113805172146428, -0.140173723617088, -0.218646365692892, 0.145251404907778, 0.241912200385107, -0.39076691014967, 0.269582742739157, -0.114902811187521, 0.0334101383852755, -0.00686857701195878, 0.00100420142513388, -0.000102561832459519, 6.97559986676052e-06, -2.84442442190961e-07, 5.26835606386175e-09, -5.2386862509022e-301, -4.23361579291131e-300, 3.51194471524504e-300, 8.48673232352492e-300, -4.71477653411186e-299, -9.34425566958953e-299, 1.75113216586752e-298, -2.92773023008038e-298, -2.47570628667144e-298, 3.39914400391618e-297, 1.37771737014351e-297, -9.89264407246447e-297, -5.91538277540067e-297, -1.77857798477124e-300, 0.0034120255431575, 0.0132147181054833, 0.030880834016168, 0.0559953701360432, 0.0844404638849241, 0.106297954454242, 0.105492275401956, 0.0648803301052844, -0.0195397564076787, -0.11794292483858, -0.158622610937647, -0.0690492778347738, 0.119332206147761, 0.200768997230937, -0.0140794793119891, -0.259789642949025, 0.0223307664223831, 0.323984679012347, -0.367844013054059, 0.217340710940809, -0.0834695297879435, 0.0223952028111374, -0.00430897062566834, 0.000595275883078205, -5.78518951802092e-05, 3.76438326621822e-06, -1.47490448466161e-07, 2.63417803193088e-09, 1.92319351518369e-300, -4.6995987830483e-301, -1.10182368766949e-299, -6.43023162649834e-300, 9.88937756544796e-301, -1.67716831232187e-298, 2.36432062318119e-298, 4.65722610220725e-299, -3.10307308552417e-297, 2.05064140382095e-297, 2.45807862155749e-296, 1.09786817924451e-296, -0.0270566480764053, -0.0468634891508447, -0.0596602195596015, -0.0666139703071929, -0.0644452220604183, -0.0479622399421162, -0.012645443962329, 0.0405297976454582, 0.0984833209915629, 0.131952125543185, 0.103504300479784, -0.00234076256049541, -0.134523725173483, -0.167232288490688, -0.00791120502378967, 0.20321252444226, 0.11862447451755, -0.232143202573217, -0.10595572976593, 0.368916223962043, -0.328096907259969, 0.169930465717401, -0.0594172393779041, 0.0148003595813558, -0.00267618459028862, 0.000350403116692284, -3.248012246049e-05, 2.02575184274583e-06, -7.63798078664119e-08, 1.31708901596544e-09, 3.56327789603186e-300, 7.20276786034171e-300, -1.71209380162471e-299, -2.36268117258936e-299, 1.63215019512819e-298, -1.97389958062255e-298, -1.28539236843526e-298, 3.1439820242178e-297, -1.13729638125531e-297, -1.99780848678712e-296, 4.54434467063051e-298, -2.75473790449479e-300, -0.0029781978737824, -0.0115345107668696, -0.0269982963171618, -0.0492238525966123, -0.0752100628045356, -0.0974028126722077, -0.103024239289136, -0.0769798936859321, -0.0112731091629928, 0.0793322333432672, 0.145452119789026, 0.118717318548918, -0.0233842051656483, -0.17487011180737, -0.133728957018368, 0.126720259702737, 0.21611868587915, -0.149650717727931, -0.218083190834032, 0.380176180348483, -0.280220866943655, 0.129389846509727, -0.041539002798914, 0.0096565016144348, -0.00164691689374471, 0.000204933724325244, -1.81567597578199e-05, 1.08729028832708e-06, -3.95071822269082e-08, 6.58544507982719e-10, -4.11660917938333e-300, 6.46902324442826e-299, 1.24750198814825e-298, -9.91316791820091e-299, 2.72686762185217e-298, 3.48815376243921e-298, -3.19453952680753e-297, -2.2766636868257e-297, 2.55019299408845e-297, -1.36625059763015e-296, 0.0253376470002428, 0.0438860919486656, 0.0559671671950235, 0.0629577752790979, 0.062269585028844, 0.04956520626586, 0.0206670030172405, -0.0248994148287453, -0.0787651593135741, -0.119515488410533, -0.115678169661124, -0.0439878985573488, 0.0769036016271966, 0.162530321218288, 0.10416764096547, -0.0935922024464643, -0.204125696241691, 0.00458947759712971, 0.254609727178844, -0.035583413655204, -0.301220357447554, 0.365072930581766, -0.230855363053603, 0.0962635781645432, -0.0285738611712296, 0.00622735749352191, -0.00100500903267829, 0.000119143394447565, -1.01091688837557e-05, 5.821694545906e-07, -2.04122241533396e-08, 3.2927225399136e-10, -2.88344259626884e-299, -1.22673122331788e-298, 3.92878221079514e-299, -3.30025081611362e-298, -2.77580850767713e-298, 2.845135713326e-297, 1.53378858871487e-297, -3.47285478457499e-297, 1.09516160100926e-296, -4.10652362160783e-299, 0.0026221914133447, 0.0101557036744521, 0.0238025822222568, 0.0435905392010923, 0.067307636998915, 0.0891251157736241, 0.0988252304934229, 0.0834206112161391, 0.0334309666893541, -0.0455081157623019, -0.120882091540368, -0.136307454345572, -0.0485398053520555, 0.105207292864011, 0.177316968310484, 0.0328617026309607, -0.195465430289109, -0.120427082859639, 0.230358552285757, 0.0857752016463125, -0.350883157930984, 0.332134924570163, -0.184441650469399, 0.0701640761635409, -0.0193702476636685, 0.003973400099006, -0.000608556064657587, 6.88856522437663e-05, -5.60751591119739e-06, 3.11005444185247e-07, -1.05354258294671e-08, 1.6463612699568e-10, 3.89104525375912e-299, -3.77645086687371e-299, 2.23482330559889e-298, 9.53726444920212e-299, -2.48489551035971e-297, 1.86355515762348e-297, 1.75308522261542e-296, -6.78980864606906e-297, -0.0238240825235075, -0.0412645213744287, -0.0526994479378465, -0.0596437482396808, -0.0600420996705081, -0.0502499562021839, -0.0264811655050389, 0.0125444694521613, 0.0615978397043627, 0.105079319015526, 0.117419763793362, 0.0738763132375494, -0.0253840765530702, -0.129562451270511, -0.144948627057307, -0.0146398830963567, 0.162197086629623, 0.147059686598357, -0.116751590365643, -0.212079838893899, 0.15518807448136, 0.194753426728665, -0.368938551280933, 0.289329473998083, -0.143498158564273, 0.0502124893718115, -0.0129580261658406, 0.00251059900805333, -0.000365863978895924, 3.96241972765515e-05, -3.09965072143924e-06, 1.65792360586762e-07, -5.4323685333692e-09, 8.23180634978399e-11, 2.24885069306622e-299, -8.49213545018114e-299, -2.10305709256661e-298, 2.44641201681067e-297, -1.2203750429108e-297, -1.38175846563452e-296, 1.88339173406816e-296, -1.78359774378221e-298, -0.00232643782873685, -0.00901025496668477, -0.0211411518868607, -0.0388584183255136, -0.0605179512029141, -0.0815697383111159, -0.0937861413435018, -0.0861132662822986, -0.048878871018039, 0.0176734976934142, 0.0930771465535938, 0.133833774376033, 0.0928757697984479, -0.0324138298293331, -0.154057220148177, -0.128810482030977, 0.071155802030853, 0.204651937593558, 0.000917019193115771, -0.249614874103932, 0.0492637905594502, 0.278836145158064, -0.361105130868911, 0.243107144405384, -0.109071866772257, 0.035346750391624, -0.00856413234964062, 0.00157211513491868, -0.000218501632684761, 2.26839141862659e-05, -1.70781453290448e-06, 8.82056392174331e-08, -2.798511502623e-09, 4.115903174892e-11, -5.93642606072246e-299, 5.23781528149188e-298, -2.36166753867305e-297, -2.20859366156033e-297, -3.50700369218867e-297, -3.25813543762601e-296, 0.0224811952175471, 0.0389385723316661, 0.0497884330594682, 0.0566337353708313, 0.0578355551606846, 0.0503138919507592, 0.0306728899129835, -0.00278968739309491, -0.0469486098780287, -0.0904779955930191, -0.11297774747345, -0.0911312065184359, -0.0153940868705545, 0.087440945666031, 0.146759245032523, 0.088037063579222, -0.0748694550445199, -0.180280640962876, -0.0518321274057499, 0.188681460743756, 0.119495353940108, -0.229403602657931, -0.0658079623280711, 0.332741395412873, -0.334697913032057, 0.198088265106876, -0.0812078561346212, 0.0245133493927182, -0.00559765643503988, 0.000976291473231472, -0.000129690704056448, 1.29282256992322e-05, -9.38087883931937e-07, 4.68398687555089e-08, -1.44041910716965e-09, 2.057951587446e-11, 5.11407664667893e-299, 1.56073136489162e-297, 1.87262443386301e-297, 3.57134492696913e-297, 3.01067835126666e-296, -3.61629804671805e-298, 0.00207806285254703, 0.00804830282028694, 0.0189014745758805, 0.0348478692787147, 0.0546583983691998, 0.074744301374297, 0.0884323005687807, 0.0863557852051555, 0.0592859037334796, 0.00439362123176255, -0.0662001073963658, -0.120294103912396, -0.114591417300028, -0.0265481567894786, 0.101993443726285, 0.157538593267011, 0.0457361676044914, -0.148628062064442, -0.155742235113205, 0.10977030469273, 0.206966956502187, -0.161362486106532, -0.171879875475703, 0.357002281444955, -0.296970494341223, 0.15716356653323, -0.0593527842665138, 0.0167706107460415, -0.00362151686052126, 0.000601624941919007, -7.65357433872048e-05, 7.33749012838913e-06, -5.1380958474164e-07, 2.48296187252992e-08, -7.40791111381866e-10, 1.028975793723e-11, -2.24066654096464e-297, 2.57409910381798e-297, 9.88217834946605e-297, -2.51843502854194e-296, -0.0212816502681797, -0.0368608995333991, -0.0471793262444005, -0.0538928417801408, -0.0556915912862106, -0.0499582735431423, -0.0336631066209816, -0.00491156907340564, 0.0345857622915974, 0.0766380627718875, 0.105131198138009, 0.0993307865079841, 0.0452185883455059, -0.0463278736087167, -0.126520898685543, -0.124501302937105, -0.00876149374982767, 0.139009840017824, 0.145430724917558, -0.0523901068506406, -0.203641230091639, -0.00334229155239343, 0.244704412764028, -0.0630482839383654, -0.256722185146756, 0.355963783728098, -0.254118863941888, 0.121807850318846, -0.0426606710224011, 0.0113315026562792, -0.00232099260135682, 0.000368093455385365, -4.49246874345975e-05, 4.14816878940934e-06, -2.80666439932429e-07, 1.31402313931041e-08, -3.80686279394583e-10, 5.14487896861499e-12, -1.67116118964703e-297, -8.38886287388434e-297, 3.66721376073976e-296, 9.71378484043074e-297, -0.00186745723208211, -0.00723263075960861, -0.0169990960812777, -0.0314211564101593, -0.0495772040034796, -0.0686104779796191, -0.0830695156534964, -0.0850243645688121, -0.0659726096787421, -0.0214395641429474, 0.0421206521870995, 0.101605177089896, 0.120244548014451, 0.0676016392702968, -0.0451760544869131, -0.141525300570246, -0.116379318976619, 0.0465345033374857, 0.179019170607708, 0.0659803259590539, -0.183093383809196, -0.116520274276821, 0.228916206813574, 0.0461698185797689, -0.314408655384955, 0.335841052092188, -0.210846353029052, 0.0924607554307679, -0.030201204484333, 0.00756926647351932, -0.0014745354282645, 0.000223710680151413, -2.6237311633637e-05, 2.33651175827045e-06, -1.52924036405232e-07, 6.94311662653947e-09, -1.95488476130819e-10, 2.5724394843075e-12, -6.8191923006481e-297, -4.99609618976687e-296, 0.0202036552227493, 0.0349937573444061, 0.0448278912869238, 0.0513899921860139, 0.0536329394347676, 0.0493210453073293, 0.0357573145785533, 0.010992059884248, -0.0242158482710665, -0.0639779763357637, -0.0956383738662662, -0.101335837201797, -0.0655768663340076, 0.0107490560845645, 0.0964823503657013, 0.132699615167439, 0.0694679340813669, -0.0701538332661074, -0.160636143131251, -0.0703429566975847, 0.135632305247937, 0.161147047998105, -0.105281856140226, -0.201029929077018, 0.167807972508857, 0.149442883182244, -0.344346547072583, 0.303190113760121, -0.170325227739901, 0.0688854684689907, -0.0210864033927562, 0.00500299350173399, -0.000929183622004629, 0.000135114242016625, -1.52511118944049e-05, 1.31152186994803e-06, -8.31227125340201e-08, 3.66320372955605e-09, -1.00316894550819e-10, 1.28621974215375e-12, 5.21422136476317e-296, -1.05321267453641e-297, 0.00168733084269139, 0.00653500425328589, 0.0153696509558603, 0.0284714545421937, 0.0451493331926984, 0.0631105739479152, 0.0778709985539559, 0.0827090663057186, 0.0699448176210601, 0.0343276257536381, -0.0214721953559171, -0.0813637555126592, -0.115535201307353, -0.0923176745562691, -0.00450055827417475, 0.103351054155752, 0.140931629940339, 0.0438202707296069, -0.121190000812142, -0.156191952971563, 0.0371393373836287, 0.201851553264775, 0.00336347342362027, -0.239763246642033, 0.0767012743188446, 0.234815863580925, -0.34967810749133, 0.263903340306895, -0.134391222081947, 0.050461913115782, -0.0145364229684276, 0.0032745657411056, -0.000581100898665069, 8.11283243032526e-05, -8.82579521794891e-06, 7.33777653934233e-07, -4.50794913787991e-08, 1.92999704583093e-09, -5.14447700924342e-11, 6.43109871076874e-13).finished();
}
} // namespace mrcpp
} // namespace detail
