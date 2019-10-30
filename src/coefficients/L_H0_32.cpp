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
auto get_L_H0_32() noexcept -> Eigen::Matrix<double, 33, 33> {
  return (Eigen::Matrix<double, 33, 33>() << 0.707106781186548, 4.91676504372299e-310, -8.75798773413159e-309, -5.53136067418837e-310, 7.25222843949142e-309, -7.3751475655845e-310, -2.538279953822e-308, 1.22919126093076e-310, -4.50498597131119e-308, 1.22919126093076e-310, 1.25070210799704e-307, -1.84378689139614e-310, -6.84229315397101e-307, -1.22919126093076e-310, 9.63931786821893e-306, 6.14595630465379e-311, -7.11533955471787e-305, -6.14595630465379e-311, 4.43169667156946e-304, 0, -1.43046604438574e-303, 6.14595630465379e-311, -2.93845933268308e-303, 0, 4.573757009816e-303, -3.0729781523269e-311, -3.06689223156719e-301, -6.14595630465379e-311, -1.31082152588085e-300, -9.21893445698069e-311, -5.46799927984166e-301, 3.0729781523269e-311, -7.03153571435811e-299, -0.612372435695794, 0.353553390593274, 6.30500693167848e-309, 2.88739908109651e-309, -6.39851669654811e-309, 1.20172651869433e-308, 2.05021659954286e-308, 1.59425146237046e-308, 3.9083840683862e-308, 1.57201438413926e-307, -1.08206642997765e-307, -1.33361489976951e-306, 5.91942379763317e-307, -2.34375203410507e-307, -8.34814883929307e-306, -2.1462528085989e-306, 6.16202553161095e-305, -8.20923795890293e-305, -3.83797085532474e-304, 5.51627447613675e-304, 1.23881957045493e-303, -2.83864110510177e-303, 2.54478097918783e-303, 1.26909041210214e-302, -3.96098872712933e-303, -6.14327633931639e-302, 2.65600658763987e-301, 3.35436834993675e-301, 1.13520474142291e-300, -3.59750808935793e-300, 4.73542628400865e-301, -6.46932901144903e-300, 6.08948855630904e-299, -4.9628597160079e-309, -0.684653196881458, 0.176776695296637, -7.49806669167757e-309, 1.53648907616343e-309, -2.34775530837773e-308, -1.35211038702382e-309, -3.34954618603629e-308, -2.27400383272189e-308, -3.03978998828174e-307, 7.47962882276361e-308, 2.5800109971306e-306, -9.194350631762e-308, 4.54432009166098e-307, 4.65353373519468e-306, 4.15663316796342e-306, -3.21314283181081e-305, 1.58971121397485e-304, 2.08961100788277e-304, -1.06822239643644e-303, -6.42120357235329e-304, 5.49700373363278e-303, -1.52917084015323e-303, -2.45758311732495e-302, 4.23828447147458e-303, 1.18964036225397e-301, -1.48278931426945e-301, -6.49570637610979e-301, -6.10538253477648e-301, 6.96654446058985e-300, -5.32411127063034e-301, 1.25278017602472e-299, -3.43418247227743e-299, 0.233853586673371, 0.405046293650491, -0.522912516583797, 0.0883883476483184, -3.5985054316584e-309, 1.72132391049753e-308, -3.76439823660042e-309, 2.17607666175828e-308, 5.48248111545292e-308, 2.19836535641938e-307, -1.83153339101372e-307, -1.7913469993795e-306, 4.46110000207328e-308, 1.36445031491676e-308, -1.0580848604616e-305, -4.18762655339382e-306, 7.15114488782782e-305, -1.24903354328019e-304, -4.715505443661e-304, 8.05300477749257e-304, 1.42633799615658e-303, -4.17224494467073e-303, 3.5515463025303e-303, 1.86760823590333e-302, -1.1024387232496e-302, -9.05279167057159e-302, 3.37187066678117e-301, 5.07854825304306e-301, 1.37248303555726e-300, -5.24120116335517e-300, 1.394056272004e-300, -8.64681208687065e-300, 7.83298689331976e-299, -1.59794863920997e-309, 0.153093108923949, 0.592927061281571, -0.350780380010057, 0.0441941738241592, -1.27835891136798e-308, -5.59282023723491e-309, 1.4750295131169e-309, -8.10037040953364e-308, -8.84403112239674e-308, 2.82959828266258e-307, 4.65064513573149e-307, -2.08716676106041e-307, -1.22089421991947e-306, 1.68302096637899e-305, 5.93170826787351e-306, -1.14420559174563e-304, 8.69163598986654e-305, 7.51099716914693e-304, -4.49031864659014e-304, -2.27177039220077e-303, 2.42265011628493e-303, -5.58886531435286e-303, -1.09221936272539e-302, 1.63656930622213e-302, 5.33581124851833e-302, -5.31993093481365e-301, -3.45129625767267e-301, -2.18395265055536e-300, 2.88605523895669e-300, -2.07241757936142e-300, 2.10100282742945e-300, -1.23537745818569e-298, -0.146575492494482, -0.253876200144874, -0.163876382526586, 0.581703452155821, -0.219863238741723, 0.0220970869120796, 3.00479644957212e-309, 3.07873998636249e-309, 3.69612050327841e-308, 1.51020550990432e-307, -2.09748524551252e-307, -8.02382444437052e-307, -2.44273914245481e-307, 2.1627812297211e-306, -1.23713092418597e-305, -1.04992466300867e-305, 7.76168399507099e-305, -1.5338418735387e-304, -5.35461767314638e-304, 7.86475440956814e-304, 1.47899396703274e-303, -4.20672559804505e-303, 4.31451297302661e-303, 1.89869979921649e-302, -1.55734749045141e-302, -9.27406324674985e-302, 3.73415290400173e-301, 6.05300758242745e-301, 1.51501636279003e-300, -5.01861076363182e-300, 2.03300142661262e-300, -3.43221781854237e-300, 8.79608283498267e-299, -4.73699582181187e-308, -0.0689981317681863, -0.26722861525761, -0.421585548851001, 0.478033079399324, -0.132121363478811, 0.0110485434560398, -4.19154219977385e-308, 2.20025235706604e-308, -2.55118646206177e-307, 1.1873987580591e-307, 1.66155928696314e-306, 1.17504538588675e-306, -2.10191705619158e-306, 6.40543857983623e-306, 1.26259453344654e-305, -2.43681021523216e-305, 2.1851959933109e-304, 2.38048006731749e-304, -1.18453308301293e-303, -3.13159029381746e-304, 6.07426077430093e-303, -2.78809654498371e-303, -2.74373252482655e-302, 1.75048543603299e-302, 1.33536559816293e-301, -1.55953280585873e-301, -8.51943752256502e-301, -5.741135526019e-301, 7.43380786337768e-300, -2.39295254001037e-300, 6.84404222606757e-300, -3.99064155879635e-299, 0.106977062012728, 0.185289706650491, 0.181798066847189, -0.0566069404148025, -0.534885310063639, 0.3548027758708, -0.0771422564770762, 0.0055242717280199, 1.6734286650765e-308, 1.39097395759403e-307, -1.60391213743683e-307, -9.06549681661225e-307, -1.53047180081891e-306, 1.85147029708261e-306, -1.03837084133932e-305, -1.07883783429107e-305, 4.38320490347072e-305, -1.7427823920955e-304, -4.00052932532258e-304, 8.59149622661798e-304, 6.01321539781313e-304, -3.95671699760532e-303, 4.42442741162871e-303, 1.8134189094044e-302, -2.55475105391448e-302, -8.80124509270195e-302, 2.57498681476712e-301, 6.32287857262977e-301, 9.92039693896816e-301, -4.93447989359874e-300, 3.51727762613923e-300, -1.67611473276304e-300, 6.54735767160818e-299, 6.31804308118405e-308, 0.0394511911654769, 0.152793806371937, 0.301313449619951, 0.204994402553155, -0.528802957971469, 0.246372609862836, -0.0441077726196727, 0.00276213586400995, -1.63666816392929e-307, 4.21305304684014e-307, -6.60383004935045e-307, 1.15593146177928e-306, -2.39667712056277e-306, 1.57507953580036e-305, 1.08603964668275e-305, -8.18068576652285e-305, 1.2661241561523e-304, 6.52573215339595e-304, -4.2357998457193e-304, -1.25817171891102e-303, 1.23637342546837e-303, -6.28217189991412e-303, -6.12815060880528e-303, 2.83493678322793e-302, 2.93512014243342e-302, -4.07384946770096e-301, -3.80398308615442e-301, -1.70894422386013e-300, 1.59593171243243e-300, -3.976299326496e-300, -6.89238880048125e-300, -1.01662326451977e-298, -0.0842790976968415, -0.145975679226991, -0.161531821313557, -0.0637090094933597, 0.216717679791878, 0.39931734961447, -0.455808912293701, 0.163205770906627, -0.0248208301311305, 0.00138106793200498, -3.7905185508952e-307, 5.4181502384453e-307, -8.94487282107669e-307, 2.78643446329046e-306, -1.10717338171141e-305, -1.38938455567531e-305, 4.94996626792527e-305, -1.84013614533588e-304, -4.56103746732885e-304, 6.821870707751e-304, 4.48705625774695e-304, -2.18134430040774e-303, 4.93959210302965e-303, 1.02156997871292e-302, -2.17998317212625e-302, -4.80943694976608e-302, 2.23441789980385e-301, 5.54132406107742e-301, 1.10155594383179e-300, -2.82675644078605e-300, 3.33914442824588e-300, 7.17927241358248e-300, 6.1694234487978e-299, -4.18877651943676e-307, -0.0255777360424, -0.0990621457259044, -0.210503024770429, -0.253153933531557, 0.0104952076465843, 0.490608224764219, -0.359502355151361, 0.10437804074912, -0.0137934051677984, 0.000690533966002488, 2.47866417766686e-307, 1.43606415014539e-306, -1.79640156828724e-306, 4.61303188314701e-306, 1.26371309749399e-305, 3.41555375674827e-306, 2.40031429750387e-304, 1.54971271574853e-304, -1.06169926279337e-303, 8.15496616857914e-304, 3.60998477201978e-303, -3.23064989564243e-303, -1.59895393328273e-302, 1.81285460675857e-302, 7.35296532810553e-302, 3.52877535354799e-302, -7.16715148011513e-301, -1.6506359771073e-301, 4.8297362476408e-300, -3.15816511742217e-300, -2.66432168754036e-300, -6.10506941041798e-300, 0.0695453758035306, 0.120456124323186, 0.1406978842158, 0.0963808314913037, -0.0645778489604212, -0.291065800584817, -0.208957848692535, 0.493804817536819, -0.266261292339, 0.0649588425253176, -0.00758802259175784, 0.000345266983001244, 2.01458685832514e-307, -1.24594283308362e-306, -6.95224719317829e-306, -6.21750291848529e-306, 1.43912464432398e-305, -1.85288056358117e-304, -2.95462333028651e-304, 7.5289352661797e-304, -5.93079477710245e-304, -1.2119596761462e-303, 4.52404515991594e-303, 4.33605203771478e-303, -2.10810440493576e-302, -1.31209144619447e-302, 2.93714329070491e-302, 4.01380712691585e-301, 5.77698619640815e-301, -2.13280036944693e-300, 3.79956171860076e-300, 4.9827795684801e-300, 2.49731387453393e-299, 9.18974116453351e-307, 0.0179405987020252, 0.0694836399939341, 0.153009670967929, 0.220107701663293, 0.143140129389143, -0.180507282585682, -0.371745290494229, 0.441307394385264, -0.18812298593045, 0.0395553021143058, -0.00413960570259111, 0.000172633491500622, 2.68516830950322e-307, 1.24238662911684e-305, 3.10690383112856e-306, -6.10543601473716e-305, 1.2759103623762e-304, 5.62723943632786e-304, -2.910470463293e-304, -2.82421442493971e-304, -1.83740443869287e-303, -6.005137695419e-303, 1.04771228928928e-302, 1.66425423348912e-302, -6.32239532784928e-302, -1.8304004199228e-301, -5.09272965133244e-302, -1.46647200908946e-300, -1.61926841528603e-300, -3.4771241162372e-300, -1.1636732213173e-299, -6.30434725202139e-299, -0.0592039757166833, -0.102544293951369, -0.123357924633427, -0.103239338518338, -0.00642871453489237, 0.168093812205297, 0.274914431263867, 0.00810642466595898, -0.454921548126905, 0.363635843788664, -0.128117052855236, 0.0236610244206259, -0.00224257483775315, 8.6316745750311e-05, -7.85640283279772e-306, -4.24921047602496e-306, 4.12223789968783e-305, -1.73700504912327e-304, -4.12000063715014e-304, 4.8240064011413e-304, -1.4795320599118e-304, 9.4517956760784e-304, 4.28182773546681e-303, -8.30780948152857e-303, -8.6339588692002e-304, 5.81637391014992e-302, 2.42326669283394e-302, 1.52853280080605e-301, 1.14202778956438e-300, 7.78252603982521e-301, 1.81952071402203e-300, 1.06694695297326e-299, 2.95669725739662e-299, 7.66984617039264e-306, -0.013284298602997, -0.0514498672554539, -0.115603486198934, -0.181283661577571, -0.176520799885996, -0.00586588815455456, 0.273642327022754, 0.199320449302147, -0.467026783576276, 0.282240662668031, -0.0847113439910263, 0.0139448970444285, -0.00120766350936337, 4.31583728751555e-05, 7.2312092689295e-306, 7.11376004394757e-306, 2.15707517023896e-304, 1.30510181103642e-304, -8.23345617654587e-304, 1.18554722726276e-303, 9.69468986154571e-304, -2.08566929880678e-303, 2.4452733199828e-303, -1.01009902441289e-302, -3.71976290320219e-302, 2.27943030527659e-301, -2.74457685599638e-301, -3.60136037617216e-301, 1.28145557443796e-300, -8.87562955171961e-301, -3.33997276915254e-300, 2.37508071107309e-299, 0.0515434170219584, 0.0892758170777425, 0.109344092753528, 0.101404123890135, 0.0407301127516175, -0.0856743299651528, -0.220498364544324, -0.184736218894801, 0.159017677884524, 0.344078934089262, -0.430008059376046, 0.209163768146521, -0.0546672604778347, 0.00811598337070386, -0.000647015840060858, 2.15791864375777e-05, 1.20837644341596e-305, -1.63914473464059e-304, -2.36712459696966e-304, 5.00561498293467e-304, -8.5752342855184e-304, 7.93452020375125e-304, 2.75875527135307e-303, -1.32051796848463e-302, 1.07989071035754e-302, 9.84924695318794e-302, -1.62075073904864e-301, -7.44528369435192e-302, 7.90023171454265e-301, 1.15500793349816e-300, 1.20771063679377e-300, 1.80062280017371e-300, -4.54262161158189e-300, -6.16596600189262e-305, 0.0102345192513237, 0.0396381226168159, 0.090193190712834, 0.148631623034035, 0.172042866302963, 0.0907318468848569, -0.119866611828528, -0.278556309580142, -0.0124278691962095, 0.424348123789452, -0.366809959680973, 0.149375400964147, -0.034568341247478, 0.00467292785582249, -0.000345098354428059, 1.07895932187889e-05, 1.37366179740172e-304, 4.68960373815106e-304, 1.36971240588035e-305, -2.76506143929431e-304, -3.42484613538309e-303, -3.67090256901806e-303, 2.74671059712162e-302, -1.79067482766213e-302, -1.76524254580051e-301, -2.4470428268463e-302, 4.49290364011405e-301, -1.75430323208956e-300, -4.90772413896718e-300, -6.30798402108148e-301, -6.27125089075093e-300, -3.98420729816288e-299, -0.0456399474313251, -0.0790507078058277, -0.0979718639481214, -0.0966015606026778, -0.057448884878591, 0.0338731370251748, 0.156041624778151, 0.210533095503393, 0.0531636663628142, -0.258761013058376, -0.183866825334376, 0.444497812568922, -0.295414932464577, 0.103484562945811, -0.021484146740091, 0.0026655154492612, -0.000183343732310454, 5.39479660939444e-06, -3.97080418678482e-304, 1.36769685871635e-304, 1.57600475190368e-304, 2.08951595298608e-303, 1.60634794317233e-303, -2.39695603808034e-302, 3.677528816347e-302, 1.67142750359878e-301, -8.34202174665686e-302, -3.3135367788384e-301, 1.56697586339225e-300, 4.14358116171676e-300, -1.3293969582876e-300, 4.96780434374525e-300, 1.60633629131474e-299, 4.94134979083506e-304, -0.00812779240754105, -0.0314788046358376, -0.072235100975407, -0.122860358428221, -0.156506062959549, -0.124232982942462, 0.0142344139090215, 0.201327606218267, 0.205689264120524, -0.147671631501979, -0.317271004569918, 0.420035413429816, -0.227379255429147, 0.0698910274057256, -0.0131549524088687, 0.00150807484138255, -9.70688678741415e-05, 2.69739830469722e-06, -4.4952120569999e-304, 7.72381442729944e-304, 6.33408955849987e-304, 6.75524674563694e-304, 1.49233270943194e-302, -4.79739796078428e-302, -1.30958808318538e-301, 3.07115740816215e-301, 1.52149554866555e-301, -7.7582898723407e-301, -1.79240334815711e-300, 2.19070870391167e-300, 3.24779190918296e-300, 3.29491064655878e-299, 0.0409507954903469, 0.0709288583996433, 0.0886307272903126, 0.0909639704527409, 0.065258486931141, -0.00145260345997924, -0.103223850969413, -0.186198844713189, -0.144326934607074, 0.080035135058834, 0.275687717919289, 0.00888656633110127, -0.396916132760715, 0.368867993548947, -0.168683485161818, 0.0461917552579768, -0.00795115139099776, 0.000847091515042966, -5.1232818674329e-05, 1.34869915234861e-06, -2.53630201222814e-304, 1.41278941361967e-304, 1.72235977302046e-304, -2.19439657978358e-302, 4.34639990041421e-302, 1.73772720092888e-301, -2.27171868711545e-301, -4.39319448032024e-301, 1.01914380627417e-300, 3.74747336701846e-300, -1.21781705959992e-300, -6.51452374300709e-300, -1.35613111710366e-299, -4.99874165270136e-304, 0.00661134635589627, 0.0256056343323954, 0.0591088434232506, 0.102723890512604, 0.138970988827876, 0.133391498453093, 0.0478956046615338, -0.110238279264893, -0.219985482132738, -0.0834232810706853, 0.247772805434256, 0.165701105440151, -0.424277917160628, 0.306332923264078, -0.121360459018218, 0.0299640781761113, -0.00475142615106978, 0.000472770516148696, -2.69655523597692e-05, 6.74349576174305e-07, -1.37804220483874e-303, -2.03309746463197e-303, 3.36526831050813e-302, -4.30107511959162e-302, -2.36067977940399e-301, 2.43412427867812e-303, 7.58497020990783e-301, -1.82965990166813e-300, -7.28775128621149e-300, 9.67117534545153e-301, 2.2599098722235e-300, -3.51355356079941e-299, -0.0371360054612984, -0.0643214482491242, -0.0808534054499439, -0.085324656829441, -0.068359789048132, -0.0189633117232361, 0.0635837506744712, 0.149958365846225, 0.168394422018244, 0.0429448711120381, -0.178415943868476, -0.215433292864517, 0.142996615456625, 0.291388108825173, -0.410468654023286, 0.243326247690153, -0.0850729385293069, 0.0191236468801429, -0.00281084222623789, 0.000262346470332942, -1.41573265450511e-05, 3.37174788087152e-07, 9.07498457903922e-304, -2.8284071179698e-302, 5.69262434733852e-302, 2.04874177305757e-301, -5.25626825049685e-302, -5.60649713035165e-301, 1.59804227387892e-300, 6.429733021329e-300, -2.24235770920799e-300, -1.80416005235367e-300, 1.88776227270097e-299, 3.83495198655584e-303, -0.00548335795155571, -0.0212369540276693, -0.0492384817480257, -0.0868997867372054, -0.122464688632806, -0.131322704818521, -0.0816110359705669, 0.0384968769960545, 0.172066290983573, 0.175488770814139, -0.0501366381202787, -0.270683109786265, -0.000803764424668383, 0.371470721264963, -0.369741542142612, 0.186340754749332, -0.0583166494034511, 0.0120315623487937, -0.00164793638714457, 0.000144828248486399, -7.41592932462492e-06, 1.68587394043576e-07, 7.51817927222498e-303, -5.26193252698315e-302, -1.47364234680971e-301, 2.48319661270821e-301, 2.65897171098419e-301, -6.99050476231603e-301, -3.64099947178029e-300, 2.27070834891808e-300, 8.70715035593979e-300, 2.76907742959666e-299, 0.0339717517016025, 0.05884079996929, 0.0742936264245927, 0.0800037948889859, 0.0689062655971775, 0.0318660438635904, -0.0346253730907701, -0.114902269304468, -0.161615833327441, -0.107864276927582, 0.0640047727487351, 0.21795005162665, 0.102133119058666, -0.240190695849873, -0.146237150827973, 0.405311157739788, -0.315303356189628, 0.138384156152882, -0.0392045609971083, 0.00747405747584578, -0.000958375545155694, 7.95788306785343e-05, -3.87659371888854e-06, 8.42936970217881e-08, 4.51998313356605e-302, 1.57901421271156e-301, -1.57050665522888e-301, -4.4372949208564e-301, 6.92252394132427e-301, 4.95983213539395e-300, -8.60742758993214e-302, -1.2004276727075e-299, -1.06463754867098e-299, -7.86039375561881e-303, 0.00462151457774989, 0.0178990489938801, 0.0416365256408676, 0.0743293390968425, 0.10785776491811, 0.124342504825598, 0.0981880514689788, 0.0115281199739022, -0.113108577862976, -0.184406895713343, -0.0846367449913654, 0.157465262190036, 0.218886254237158, -0.142694154981089, -0.266312347887875, 0.400808174393301, -0.257327618379487, 0.100105534682369, -0.0259085163210665, 0.00459051433800809, -0.00055330348694272, 4.35403668581781e-05, -2.02260965120282e-06, 4.2146848510894e-08, -2.09041366670318e-301, -1.04906146313838e-301, 7.13342785085615e-301, -1.35811261012238e-300, -8.05442354792808e-300, -1.44740453227851e-300, 7.12080155232894e-300, -3.91400698700552e-299, -0.0313046309536483, -0.0542212113239122, -0.0686949483603631, -0.0751077216437837, -0.0680535455514093, -0.0399825735066795, 0.0136245225316148, 0.0848280036848168, 0.142642688279874, 0.134423434769748, 0.0213966755957234, -0.146829028303034, -0.192413047610126, 0.0287491387003486, 0.265192777737401, -0.00996953747334031, -0.347307623368017, 0.36940243095303, -0.202524168778998, 0.0707847575547383, -0.0168635954707225, 0.00279083288376033, -0.000317334210055171, 2.37299576748381e-05, -1.05346045745216e-06, 2.1073424255447e-08, 8.88443142096716e-302, -4.31654346816466e-301, 1.13126091840305e-300, 6.91269776113123e-300, 1.04364926459216e-300, -4.43809021829187e-300, 2.77825402947118e-299, -1.71596859751507e-301, -0.00394815347735416, -0.0152911326660635, -0.0356609978909838, -0.0642215183683974, -0.095241307869668, -0.115508986443972, -0.104725569006555, -0.0443521460195627, 0.0608399177151585, 0.156823522657345, 0.145874985959225, -0.0241723579915451, -0.211236840069211, -0.113052959575701, 0.235223105585243, 0.12618987380636, -0.386981898851779, 0.322516108630937, -0.154574067130241, 0.0490623442773783, -0.0108283126869318, 0.001681118939966, -0.00018090493786476, 1.2887014923139e-05, -5.47807706733102e-07, 1.05367121277235e-08, -6.01081446941935e-302, -1.57738230178669e-301, -3.67172405250314e-300, -1.55214275530173e-300, 9.50784462390574e-300, 1.79629075242517e-299, 0.0290259753654005, 0.0502744640721163, 0.0638655891657856, 0.070651871387702, 0.0664442436081886, 0.0449949520307467, 0.00159792382538831, -0.0602296998266835, -0.120294637837257, -0.138964796248224, -0.0747403231416156, 0.067365852973061, 0.183893717240485, 0.113805172146428, -0.140173723617088, -0.218646365692892, 0.145251404907778, 0.241912200385107, -0.39076691014967, 0.269582742739157, -0.114902811187521, 0.0334101383852755, -0.00686857701195878, 0.00100420142513388, -0.000102561832459519, 6.97559986676052e-06, -2.84442442190961e-07, 5.26835606386175e-09, 3.06480866818991e-301, 4.15920516763064e-300, 4.49235066051211e-300, -1.27333828279715e-299, -4.89453451498632e-300, -5.6322889846039e-301, 0.0034120255431575, 0.0132147181054833, 0.030880834016168, 0.0559953701360432, 0.0844404638849241, 0.106297954454242, 0.105492275401956, 0.0648803301052844, -0.0195397564076787, -0.11794292483858, -0.158622610937647, -0.0690492778347738, 0.119332206147761, 0.200768997230937, -0.0140794793119891, -0.259789642949025, 0.0223307664223831, 0.323984679012347, -0.367844013054059, 0.217340710940809, -0.0834695297879435, 0.0223952028111374, -0.00430897062566834, 0.000595275883078205, -5.78518951802092e-05, 3.76438326621822e-06, -1.47490448466161e-07, 2.63417803193088e-09, -6.39937019797643e-300, -6.94889922989866e-300, 7.93880745463919e-300, -4.39638706392527e-299, -0.0270566480764053, -0.0468634891508447, -0.0596602195596015, -0.0666139703071929, -0.0644452220604183, -0.0479622399421162, -0.012645443962329, 0.0405297976454582, 0.0984833209915629, 0.131952125543185, 0.103504300479784, -0.00234076256049541, -0.134523725173483, -0.167232288490688, -0.00791120502378967, 0.20321252444226, 0.11862447451755, -0.232143202573217, -0.10595572976593, 0.368916223962043, -0.328096907259969, 0.169930465717401, -0.0594172393779041, 0.0148003595813558, -0.00267618459028862, 0.000350403116692284, -3.248012246049e-05, 2.02575184274583e-06, -7.63798078664119e-08, 1.31708901596544e-09, 9.70416965247145e-300, -6.60335744051658e-300, 3.75932708565414e-299, 7.47537292736164e-300, -0.0029781978737824, -0.0115345107668696, -0.0269982963171618, -0.0492238525966123, -0.0752100628045356, -0.0974028126722077, -0.103024239289136, -0.0769798936859321, -0.0112731091629928, 0.0793322333432672, 0.145452119789026, 0.118717318548918, -0.0233842051656483, -0.17487011180737, -0.133728957018368, 0.126720259702737, 0.21611868587915, -0.149650717727931, -0.218083190834032, 0.380176180348483, -0.280220866943655, 0.129389846509727, -0.041539002798914, 0.0096565016144348, -0.00164691689374471, 0.000204933724325244, -1.81567597578199e-05, 1.08729028832708e-06, -3.95071822269082e-08, 6.58544507982719e-10, 4.33345906732346e-300, 1.12197584007898e-299, 0.0253376470002428, 0.0438860919486656, 0.0559671671950235, 0.0629577752790979, 0.062269585028844, 0.04956520626586, 0.0206670030172405, -0.0248994148287453, -0.0787651593135741, -0.119515488410533, -0.115678169661124, -0.0439878985573488, 0.0769036016271966, 0.162530321218288, 0.10416764096547, -0.0935922024464643, -0.204125696241691, 0.00458947759712971, 0.254609727178844, -0.035583413655204, -0.301220357447554, 0.365072930581766, -0.230855363053603, 0.0962635781645432, -0.0285738611712296, 0.00622735749352191, -0.00100500903267829, 0.000119143394447565, -1.01091688837557e-05, 5.821694545906e-07, -2.04122241533396e-08, 3.2927225399136e-10, -4.27449831704564e-300, -7.11378259031299e-299, 0.0026221914133447, 0.0101557036744521, 0.0238025822222568, 0.0435905392010923, 0.067307636998915, 0.0891251157736241, 0.0988252304934229, 0.0834206112161391, 0.0334309666893541, -0.0455081157623019, -0.120882091540368, -0.136307454345572, -0.0485398053520555, 0.105207292864011, 0.177316968310484, 0.0328617026309607, -0.195465430289109, -0.120427082859639, 0.230358552285757, 0.0857752016463125, -0.350883157930984, 0.332134924570163, -0.184441650469399, 0.0701640761635409, -0.0193702476636685, 0.003973400099006, -0.000608556064657587, 6.88856522437663e-05, -5.60751591119739e-06, 3.11005444185247e-07, -1.05354258294671e-08, 1.6463612699568e-10).finished();
}
} // namespace mrcpp
} // namespace detail
