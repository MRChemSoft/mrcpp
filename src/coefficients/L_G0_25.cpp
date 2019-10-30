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
auto get_L_G0_25() noexcept -> Eigen::Matrix<double, 26, 26> {
  return (Eigen::Matrix<double, 26, 26>() << 6.46454430413442e-276, -0.00394815347735416, -0.0152911326660635, -0.0356609978909838, -0.0642215183683975, -0.095241307869668, -0.115508986443972, -0.104725569006555, -0.0443521460195627, 0.0608399177151585, 0.156823522657345, 0.145874985959225, -0.0241723579915451, -0.211236840069212, -0.113052959575701, 0.235223105585243, 0.12618987380636, -0.386981898851779, 0.322516108630937, -0.154574067130241, 0.0490623442773783, -0.0108283126869319, 0.001681118939966, -0.00018090493786476, 1.2887014923139e-05, -5.47807706733102e-07, -0.0290259753654029, -0.0502744640721204, -0.0638655891657907, -0.0706518713877077, -0.066444243608194, -0.0449949520307504, -0.00159792382538844, 0.0602296998266884, 0.120294637837267, 0.138964796248236, 0.0747403231416217, -0.0673658529730664, -0.1838937172405, -0.113805172146437, 0.140173723617099, 0.21864636569291, -0.145251404907789, -0.241912200385126, 0.390766910149702, -0.269582742739179, 0.114902811187531, -0.0334101383852782, 0.00686857701195933, -0.00100420142513397, 0.000102561832459527, -6.97559986676109e-06, -3.3295583827113e-274, -0.00341202554320561, -0.0132147181056696, -0.0308808340166035, -0.0559953701368328, -0.084440463886115, -0.106297954455741, -0.105492275403445, -0.0648803301062016, 0.0195397564079512, 0.117942924840241, 0.158622610939887, 0.0690492778357557, -0.119332206149437, -0.200768997233778, 0.0140794793121703, 0.259789642952702, -0.0223307664226694, -0.323984679016971, 0.367844013059292, -0.217340710943898, 0.083469529789129, -0.0223952028114554, 0.0043089706257295, -0.000595275883086654, 5.78518951810303e-05, 0.0270566481045232, 0.0468634891995463, 0.0596602196216043, 0.0666139703764346, 0.064445222127441, 0.0479622399920791, 0.0126454439756918, -0.0405297976872679, -0.0984833210936381, -0.13195212568036, -0.103504300588019, 0.00234076256172825, 0.134523725312551, 0.167232288665694, 0.00791120503476121, -0.203212524653431, -0.118624474645869, 0.232143202814597, 0.105955729885373, -0.368916224357738, 0.328096907609571, -0.169930465897987, 0.0594172394409581, -0.0148003595970489, 0.00267618459312478, -0.000350403117063514, -5.55216351527398e-273, 0.00297819799457995, 0.0115345112347165, 0.0269982974123066, 0.0492238545937752, 0.0752100658577543, 0.0974028166310627, 0.103024243487139, 0.0769798968441298, 0.0112731096701585, -0.0793322365193674, -0.145452125700864, -0.118717323467055, 0.0233842059716259, 0.174870118899038, 0.133728962699722, -0.126720264664962, -0.216118695060072, 0.149650723559443, 0.218083200635393, -0.380176196777553, 0.280220878934315, -0.12938985202131, 0.0415390045638775, -0.00965650202410352, 0.00164691696354611, 0.0253376711917863, 0.0438861338496482, 0.055967220637311, 0.0629578354287024, 0.0622696446137336, 0.0495652539089104, 0.0206670233565535, -0.02489943770314, -0.0787652335877515, -0.119515602236457, -0.115678281390532, -0.0439879437313709, 0.0769036716605674, 0.162530476771339, 0.104167747068165, -0.0935922850769304, -0.204125897451146, 0.00458946716889139, 0.254609980793777, -0.0355834256867317, -0.301220687645099, 0.365073315736917, -0.230855603440348, 0.0962636777454037, -0.0285738906167189, 0.00622736389551831, 4.07615340968166e-272, 0.00262222898273987, 0.0101558491800939, 0.0238029233097857, 0.0435911641921833, 0.0673086032936944, 0.0891263987328078, 0.0988266609361853, 0.083421834445552, 0.0334314878549772, -0.0455087239003372, -0.120883810051155, -0.136309463636439, -0.0485406249512822, 0.105208710723119, 0.177319613189339, 0.032862435175862, -0.195468232261013, -0.120429266645516, 0.230361904563091, 0.0857772222231123, -0.35088929451609, 0.332140507282542, -0.184444700381124, 0.0701652262141781, -0.0193705634736373, 0.0238275847575564, 0.0412705874217414, 0.0527071969463618, 0.0596525278309905, 0.0600509652599395, 0.0502574388553907, 0.0264852435454876, -0.0125460275103187, -0.061606565503344, -0.105094569895564, -0.117437245307056, -0.0738880228447434, 0.0253865667174401, 0.129580855520448, 0.1449711364254, 0.0146447594891173, -0.16222003391472, -0.147085259376227, 0.11676525749162, 0.212117456706541, -0.155206828658334, -0.194797267271033, 0.369009096723284, -0.289382510243353, 0.143523947871964, -0.0502214124731064, -1.41963385926271e-271, 0.00232880492313613, 0.00901942268387206, 0.0211626692044892, 0.038898009050985, 0.0605797582311716, 0.0816534546905125, 0.0938833321219283, 0.0862043899052808, 0.0489341919343682, -0.0176848259840425, -0.09316757131567, -0.133973446031556, -0.0929846182332852, 0.0324286980057213, 0.154212020015753, 0.128970132755648, -0.0711964389503593, -0.204886415747892, -0.000982154629358816, 0.249915613595208, -0.0492225997272328, -0.279305307734376, 0.361639118224668, -0.243450602586944, 0.10922238690278, -0.0225961236072843, -0.0391376341419229, -0.0500430795806484, -0.0569239567553193, -0.058133559968081, -0.0505768593235634, -0.0308410028674675, 0.00278589575169034, 0.0471661647909464, 0.0909227209254974, 0.113558685460811, 0.0916387224018386, 0.0155513115906468, -0.0878155503626423, -0.147523753546387, -0.0886328077651188, 0.0750861133212476, 0.181268191124934, 0.052430764655716, -0.189579673047332, -0.120664716926028, 0.230579347250543, 0.067158409152741, -0.335840849469693, 0.337472095835988, -0.19965024078931, 2.16630791690787e-271, 0.00211427631100764, 0.00818855694181342, 0.0192310471834585, 0.035456617155396, 0.055617311329421, 0.0760669091970067, 0.0900231186369763, 0.0879617592812306, 0.0604874678212437, 0.004690335209427, -0.0671628177706811, -0.122383327530289, -0.116915035826905, -0.0276171637305832, 0.103336348283071, 0.160678995118546, 0.0477456198782296, -0.150669369171589, -0.160085042116989, 0.110054075251681, 0.213272071614968, -0.162213758112104, -0.181813167175177, 0.370824517025669, -0.307242895158432, -0.0222753379441317, -0.0385820170750029, -0.0493840778925314, -0.0564201265466493, -0.0583285420936712, -0.0523814906080797, -0.0354152429508527, -0.00543670207421568, 0.035810579129202, 0.0798501187295951, 0.109941099009561, 0.104443020852097, 0.0485113577772002, -0.0470106782921406, -0.131784971189076, -0.131757236181329, -0.0123770706032306, 0.143752985642657, 0.155625998994963, -0.0490474975660912, -0.21547144967862, -0.0140378186088472, 0.261215001161644, -0.0489709167452709, -0.299883050617018, 0.400486518504055, 7.59138266381557e-272, -0.00203000317426484, -0.00786216848667592, -0.0184802487339753, -0.0341682305548547, -0.053945449764227, -0.0747492680270736, -0.0907178078244327, -0.0932906310727507, -0.0732092598593114, -0.0253995646951899, 0.0436259253499411, 0.109451319546991, 0.132456056479333, 0.0785589037945361, -0.0430998845397243, -0.152674120890202, -0.134745391539872, 0.0385527093624069, 0.195939987837733, 0.0927557596525619, -0.190663491423197, -0.161945837696768, 0.24028466070497, 0.129829804043791, -0.447675942673289, -0.0233065182000768, -0.0403680736700618, -0.0517226967021603, -0.0593425959502273, -0.0620719922166026, -0.0573986195060832, -0.0422617338129916, -0.0143798012281597, 0.0255867628300361, 0.0713291169211166, 0.108951147510649, 0.118400253106448, 0.0814029901876093, -0.00327864423066677, -0.103720603217951, -0.1549349498572, -0.0956432889407972, 0.0612602714501901, 0.18511394978551, 0.113006196868798, -0.127566342505531, -0.220161339214708, 0.0579377731011112, 0.290029933107368, -0.0824138518533415, -0.436474705664632, -1.54291556615051e-271, 0.00201370157784339, 0.00779903267521904, 0.0183478950647245, 0.0340212320569072, 0.0540690599145496, 0.0759084087619624, 0.0944260046720905, 0.101852438994295, 0.0890688487073363, 0.049253255595087, -0.015991175227661, -0.0900354539201532, -0.139550017461692, -0.125572238903881, -0.0314521974310528, 0.103485546040637, 0.179209515068017, 0.100096481132974, -0.102733339307652, -0.219913047872591, -0.0539818840138348, 0.236507480051358, 0.180236635754828, -0.263976638809099, -0.377800169290526, 0.0249018599800082, 0.0431312866883403, 0.0553222315978258, 0.0637542211516378, 0.0674969190085799, 0.0642644634785557, 0.0511189986559192, 0.025538332440693, -0.0128189455011132, -0.0597592390694315, -0.104435089168304, -0.129287447919859, -0.114620458878092, -0.0495852183775327, 0.052672782373495, 0.145545699935024, 0.161855060923484, 0.0604662826572541, -0.111990325100796, -0.209286895867403, -0.0957851354098003, 0.161236777060296, 0.254748081166217, -0.0214038050053274, -0.366175224470659, -0.294489766193492, -2.21932201208898e-271, -0.00198342718209194, -0.00768178044465719, -0.0180887816078185, -0.0336425600834557, -0.0538376911712282, -0.0766133395836298, -0.0977069303030655, -0.110360340417158, -0.106014500532115, -0.0767507898349184, -0.0198965904410928, 0.0560187513436661, 0.127050439462346, 0.156939088694327, 0.113131926817536, -0.00513756839271866, -0.142419718905688, -0.197456060428056, -0.0930693285337186, 0.1227245375841, 0.252621487800603, 0.111628895234481, -0.215899102254289, -0.379947250378203, -0.208510869292146, -0.0267741090926675, -0.0463741172758919, -0.0595487014969085, -0.0689443320976489, -0.0739087134895619, -0.072452619141302, -0.0618569351032695, -0.0394812335109348, -0.00410627133441806, 0.0421689549000314, 0.091905670543394, 0.131292198558983, 0.141908232330705, 0.107153949674016, 0.0235759473675818, -0.0868521304088828, -0.173233591218458, -0.173561066047275, -0.0579209496479486, 0.123682885753755, 0.238752318103653, 0.159082522839499, -0.0989396309688409, -0.330815045426613, -0.327278380013624, -0.134549498436092, 1.41602647492823e-271, -0.00192553026817606, -0.00745754666126418, -0.017577848504724, -0.0327961269894397, -0.0528609670319306, -0.0762746902744936, -0.0997351822893007, -0.117750156231692, -0.122843377964177, -0.106908134846466, -0.0641832568422879, 0.00422657159286474, 0.0856001123265791, 0.153689219748224, 0.174008813407171, 0.119470676790436, -0.00646919734940594, -0.152445212663763, -0.22913873337166, -0.160326352917365, 0.045362004016745, 0.266464220812563, 0.346718604504964, 0.242811950123939, 0.0790677680665215, -0.029136930222687, -0.050466643522283, -0.0648804501789536, -0.0754816342738231, -0.0819604575862673, -0.0827067156041902, -0.0753385896242481, -0.0572695367508016, -0.0266294470049095, 0.0163683492172946, 0.0677677935913122, 0.118460320195927, 0.154050251722882, 0.157353856793646, 0.114582405729986, 0.0247097668863172, -0.091481461251122, -0.189594045851228, -0.214744039037207, -0.131231361920667, 0.0440117064108037, 0.231911527879693, 0.330221298064422, 0.288766714177509, 0.157385662444184, 0.0421403769935433, -5.20030232942002e-272, -0.00183244531709043, -0.00709703019592699, -0.0167451169363972, -0.031345816743855, -0.0508994108792874, -0.0744918489856372, -0.0998620533355003, -0.123013569946427, -0.138102555351555, -0.137937005972827, -0.11545448771686, -0.0663891711715826, 0.00713931664113639, 0.0932319058059846, 0.169137250703955, 0.2051455108042, 0.175061288926889, 0.0713850486846124, -0.082045480184599, -0.231068023311428, -0.313501013043465, -0.295980852630426, -0.200037981758227, -0.0893849813129579, -0.0202170071135444, 0.0324537946805142, 0.0562116212850592, 0.0723554315345925, 0.0846018263992847, 0.0930753145445981, 0.0966523671880395, 0.0934448316931598, 0.0811621955626371, 0.0576462429316983, 0.0216875770950987, -0.0258796849032302, -0.0809756525935112, -0.135371173773592, -0.176737819204432, -0.190417681207495, -0.163416419009863, -0.0902528932364599, 0.0211769010996159, 0.14709078703203, 0.252346695410914, 0.303261310748389, 0.28448460245979, 0.209848866153311, 0.116827705780185, 0.0442425309621534, 0.00862764240587894, -1.09051325022612e-272, 0.00168265586849487, 0.00651689815607881, 0.0153926514072314, 0.0289131283359426, 0.0473094493437355, 0.0702425940092438, 0.09653032419024, 0.123852096234924, 0.148525329404642, 0.165496907740664, 0.168733500604349, 0.152187233335147, 0.111410742683124, 0.0456546193458447, -0.0400964467074072, -0.134098904454705, -0.218972686378082, -0.275718111074198, -0.290160761557504, -0.259535478309835, -0.195497071096424, -0.120408225035204, -0.0572819264711283, -0.0188152197395649, -0.00321469577166984, 0.0380696145832726, 0.0659385066827932, 0.0849849975123704, 0.099887094868037, 0.111369923265805, 0.118966274395121, 0.121516223455311, 0.11739098718887, 0.104724744272157, 0.0817467294716235, 0.0472494288343126, 0.00117866539486764, -0.0547375283429396, -0.11649128936729, -0.177680523987388, -0.230058311063816, -0.264935372008092, -0.275358361661381, -0.258558266336514, -0.217720415819155, -0.161970242923732, -0.103934116655383, -0.0554097103776788, -0.0230732353339541, -0.00668739335886203, -0.00101414792730299, 1.08613635665002e-273, 0.00140902355345903, 0.00545712475696082, 0.0129037627952355, 0.0243244297226091, 0.0401151913595738, 0.0604379700849314, 0.0851276955301884, 0.113576137835323, 0.144613700701483, 0.176422118719023, 0.206523734609954, 0.231901021434224, 0.249294743020444, 0.255700667513923, 0.249025997206888, 0.228781092023807, 0.196592571899811, 0.156277254616557, 0.113273086506898, 0.0734240768358654, 0.0414303262465286, 0.0195596198550282, 0.00725222072795965, 0.00187955477152005, 0.000255718835375909, -0.0546516892274594, -0.0946595024614243, -0.122165305255065, -0.144360578913904, -0.163158941359957, -0.179209482930195, -0.192610700984079, -0.203141491055956, -0.21037390685589, -0.213759566920906, -0.212722414544927, -0.206770334203257, -0.195625718414671, -0.179362896044864, -0.158527342224764, -0.134200327983991, -0.107968533806597, -0.0817680700989313, -0.0576015465111944, -0.0371736579163028, -0.0215414948538966, -0.0109040525848133, -0.00463167610230922, -0.00154891377109979, -0.000362803388423976, -4.46833039188023e-05).finished();
}
} // namespace mrcpp
} // namespace detail
