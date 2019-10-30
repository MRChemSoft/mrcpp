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
auto get_L_G0_36() noexcept -> Eigen::Matrix<double, 37, 37> {
  return (Eigen::Matrix<double, 37, 37>() << -0.0212816502681797, -0.0368608995333991, -0.0471793262444005, -0.0538928417801408, -0.0556915912862106, -0.0499582735431423, -0.0336631066209816, -0.00491156907340564, 0.0345857622915974, 0.0766380627718875, 0.105131198138009, 0.0993307865079841, 0.0452185883455059, -0.0463278736087167, -0.126520898685543, -0.124501302937105, -0.00876149374982767, 0.139009840017824, 0.145430724917558, -0.0523901068506406, -0.203641230091639, -0.00334229155239343, 0.244704412764028, -0.0630482839383654, -0.256722185146756, 0.355963783728098, -0.254118863941888, 0.121807850318846, -0.0426606710224011, 0.0113315026562792, -0.00232099260135682, 0.000368093455385365, -4.49246874345975e-05, 4.14816878940934e-06, -2.80666439932429e-07, 1.31402313931041e-08, -3.80686279394583e-10, -6.82931727520293e-255, 0.00186745723208211, 0.00723263075960861, 0.0169990960812777, 0.0314211564101593, 0.0495772040034796, 0.0686104779796191, 0.0830695156534964, 0.0850243645688121, 0.0659726096787421, 0.0214395641429474, -0.0421206521870995, -0.101605177089896, -0.120244548014451, -0.0676016392702968, 0.0451760544869131, 0.141525300570246, 0.116379318976619, -0.0465345033374857, -0.179019170607708, -0.0659803259590539, 0.183093383809196, 0.116520274276821, -0.228916206813574, -0.0461698185797689, 0.314408655384955, -0.335841052092188, 0.210846353029052, -0.0924607554307679, 0.030201204484333, -0.00756926647351932, 0.0014745354282645, -0.000223710680151413, 2.6237311633637e-05, -2.33651175827045e-06, 1.52924036405232e-07, -6.94311662653947e-09, 0.0202036552227493, 0.0349937573444061, 0.0448278912869238, 0.0513899921860139, 0.0536329394347676, 0.0493210453073293, 0.0357573145785533, 0.010992059884248, -0.0242158482710665, -0.0639779763357637, -0.0956383738662662, -0.101335837201797, -0.0655768663340076, 0.0107490560845645, 0.0964823503657013, 0.132699615167439, 0.0694679340813669, -0.0701538332661074, -0.160636143131251, -0.0703429566975847, 0.135632305247937, 0.161147047998106, -0.105281856140226, -0.201029929077018, 0.167807972508857, 0.149442883182244, -0.344346547072583, 0.303190113760121, -0.170325227739901, 0.0688854684689907, -0.0210864033927562, 0.00500299350173399, -0.000929183622004629, 0.000135114242016625, -1.52511118944049e-05, 1.31152186994803e-06, -8.31227125340201e-08, -8.09879798977923e-253, -0.0016873308426914, -0.0065350042532859, -0.0153696509558603, -0.0284714545421937, -0.0451493331926985, -0.0631105739479153, -0.0778709985539561, -0.0827090663057187, -0.0699448176210602, -0.0343276257536382, 0.0214721953559171, 0.0813637555126593, 0.115535201307353, 0.0923176745562693, 0.00450055827417476, -0.103351054155752, -0.14093162994034, -0.043820270729607, 0.121190000812142, 0.156191952971563, -0.0371393373836287, -0.201851553264776, -0.00336347342362028, 0.239763246642033, -0.0767012743188447, -0.234815863580925, 0.349678107491331, -0.263903340306895, 0.134391222081948, -0.0504619131157821, 0.0145364229684277, -0.0032745657411056, 0.000581100898665071, -8.11283243032527e-05, 8.82579521794893e-06, -7.33777653934235e-07, -0.0192296228181945, -0.0333066837314988, -0.042698053287058, -0.049097890699211, -0.0516707722203636, -0.0484977532771248, -0.0371792497098095, -0.0157911072977348, 0.0155457904503814, 0.0526380508924824, 0.0855808567463033, 0.0992729012188743, 0.0783841113160788, 0.0178836705360252, -0.0642234460082294, -0.122684026104365, -0.104571704843419, 0.00284395907453556, 0.126696024980399, 0.135325722681638, -0.0214883625107624, -0.175452197975535, -0.0762991960617233, 0.17871085801054, 0.111974246588569, -0.228614791268659, -0.0269517170439048, 0.295843084778112, -0.335612382524745, 0.222689558582825, -0.103840242456124, 0.036402053216526, -0.00990425563579717, 0.00212385753714278, -0.000360844405267845, 4.8445794222263e-05, -5.0861352982053e-06, -3.26660076033752e-251, 0.00153207183591127, 0.00593368870567777, 0.013963400136211, 0.0259150123113025, 0.0412720630161655, 0.0581811824196284, 0.0729296222870203, 0.079806636388654, 0.0719540357199555, 0.0438681329322749, -0.00426141187612788, -0.0616093793023306, -0.104741554069238, -0.104112030838646, -0.0427409109976152, 0.0587322652797063, 0.132712680318901, 0.101403537239503, -0.0384929813281587, -0.157613801987661, -0.0893970805958284, 0.123930799361527, 0.16422228181789, -0.10284410519225, -0.194414007831007, 0.174250886019147, 0.127441062768971, -0.33097242519704, 0.308024984652913, -0.182925254377905, 0.0787385774317082, -0.0258929108616281, 0.00667537840092512, -0.00136588223291527, 0.000222587496888361, -2.87801523864583e-05, 0.0183452031876897, 0.031774823996253, 0.0407601174350275, 0.0469927251718274, 0.0498091945102168, 0.0475551250194664, 0.0380942405609733, 0.0195739547585081, -0.00830860164139181, -0.0426114404664298, -0.0756028404338687, -0.0946602445399253, -0.0854511544606381, -0.0397070616886216, 0.0339002212441916, 0.102649346901492, 0.118299221053241, 0.0497649717361005, -0.0724253668581674, -0.145706237195186, -0.0721553147889438, 0.102908432536954, 0.16280393398579, -0.025046379336989, -0.199715299237934, -0.00149213179006464, 0.23468216944613, -0.0900446824196242, -0.213087933112456, 0.342282649979437, -0.272467529128813, 0.146746585206242, -0.0586961526299659, 0.0181810978224524, -0.0044540663874198, 0.000871481627204025, -0.000136448789798973, 6.33283238762327e-250, 0.00139730334054543, 0.00541173256753244, 0.0127414229918731, 0.0236854976155972, 0.0378609485089938, 0.0537600005857944, 0.0682896472495227, 0.0765831621939532, 0.0725561689866465, 0.0507591666812172, 0.00979466649249121, -0.0433990809494206, -0.0908319925103989, -0.106525407014304, -0.0692506278857319, 0.0166789461458005, 0.105898392407525, 0.125841368334886, 0.0348942344454637, -0.106869332375298, -0.146966957784035, 7.65998167320579e-05, 0.170810617362633, 0.0835867204344973, -0.175438531934332, -0.106188913000887, 0.228279874795677, 0.00822984585159455, -0.277031448659493, 0.334056276382961, -0.233590750095956, 0.11526789226928, -0.0430812061640558, 0.0126145946657282, -0.00294418044101227, 0.000551921797966277, 0.0175385735370934, 0.0303777004585288, 0.0389894342700327, 0.0450537881195608, 0.0480480314277396, 0.0465400433751737, 0.0386255088584859, 0.0225480520843916, -0.00227174380040294, -0.0338181165916788, -0.06607040173164, -0.0885536746807006, -0.0883022764197224, -0.0555258837303223, 0.00746609757690518, 0.0783492895965641, 0.116648738059885, 0.0843733227951505, -0.0170104242535956, -0.120305107736511, -0.122969984493113, 0.0093590512902203, 0.150958711606466, 0.103915475933251, -0.113821494228441, -0.165617999463698, 0.102078096441891, 0.187205826876241, -0.180482548132842, -0.105887049763435, 0.316897164622939, -0.31150634646349, 0.194908650810533, -0.0888420023610615, 0.0311734929738418, -0.00865603772866425, 0.00192915868019591, -6.78155934057293e-249, 0.0012795727714417, 0.00495576403405419, 0.0116729317556831, 0.0217299018218982, 0.0348464406947587, 0.0497890792728856, 0.0639662356214511, 0.0732163082697336, 0.072160401849457, 0.0555766140651474, 0.0210888357426906, -0.0271967353584279, -0.0757557447041825, -0.10257640293414, -0.0855112695914955, -0.0185689875334342, 0.0713686049108494, 0.124881912616841, 0.0851842144541999, -0.0393941749856511, -0.142330443778783, -0.0942259275115786, 0.0854748451918562, 0.166547228677726, -0.0157040871321673, -0.197467534879223, 0.00188014887611768, 0.229366140627104, -0.102939647760883, -0.191533580346218, 0.333817393486981, -0.279815830313235, 0.158802941767829, -0.0673085106930183, 0.0222634034096236, -0.00587887185442609, -0.0168001232416845, -0.0290986670280163, -0.0373658644859998, -0.0432636430927296, -0.0463851712766815, -0.045486154126001, -0.0388661403278003, -0.0248763961944078, -0.00276253966827742, 0.0261464634829848, 0.0571759301093, 0.0816734143986777, 0.0881480939385819, 0.0663526848929312, 0.0144278102295042, -0.0534460641611751, -0.105158877388886, -0.102565968991171, -0.029564327982924, 0.0778839357015416, 0.133785697533994, 0.0661402763840758, -0.0851740373079984, -0.153222333691268, -0.0178253184064177, 0.165862979792597, 0.0884665715548241, -0.173138918789029, -0.0994062344035534, 0.227743343439729, -0.00992811816859875, -0.257986418065618, 0.331221661799588, -0.243527457225247, 0.126670798944672, -0.0502021699122944, 0.0157086490800742, -4.27217127049867e-248, -0.00117623469312554, -0.0045555373777066, -0.0107342991881124, -0.0200073144471283, -0.0321739920311494, -0.0462202930521158, -0.0599629069700517, -0.0698301092143494, -0.0710737612504382, -0.0587898035137632, -0.030035005778189, 0.0131236991327347, 0.0607350715694515, 0.094618171651282, 0.0935792613372414, 0.0456900439738646, -0.0361979168492521, -0.107676593541995, -0.110951811581679, -0.0223574259535877, 0.100516017297782, 0.135861360325412, 0.0162329842741345, -0.142417204027014, -0.114837814732744, 0.105378340053814, 0.16579811774145, -0.102670795900445, -0.179484726886229, 0.186358790994006, 0.0848248033451671, -0.302196834430344, 0.313707164148583, -0.206251472563337, 0.0991403266047855, -0.0369125927472437, 0.0161300321044109, 0.0279380351325567, 0.0358907586597953, 0.0416283773562868, 0.044839924866299, 0.0444401376559418, 0.0389068614138357, 0.0267015744042902, 0.00696471285419338, -0.0194840913832924, -0.0490276565550923, -0.0745376208877569, -0.0859604097071582, -0.0732231240374708, -0.031862075796136, 0.0300741733172846, 0.0882232435335127, 0.107820444383065, 0.0636771387838362, -0.0319871604149683, -0.11671405670903, -0.109953720214356, 0.00770296708341166, 0.133786525259662, 0.110997042315112, -0.0695842991798577, -0.168459919693265, 0.00869370653451208, 0.19534877597099, -0.00640498826372727, -0.223904762898549, 0.115300105198725, 0.170374752716608, -0.324633324029673, 0.286201325941804, -0.170636072384865, 0.0763059068545835, 1.60595075050804e-247, -0.00108735891379097, -0.00421132296446263, -0.00992654647775794, -0.018521971761651, -0.0298588172677459, -0.0430982360918472, -0.056387758879096, -0.0666441750617077, -0.0696683053800548, -0.0609038990025437, -0.0371092906849409, 0.00109916794027439, 0.0465697067722465, 0.0845440180803185, 0.0956922783727849, 0.0651031899104497, -0.00435199216550215, -0.0821156990306545, -0.116698601640348, -0.0682708401782464, 0.0451624540316293, 0.13196414833961, 0.0909179047654874, -0.063641082020206, -0.156088988574846, -0.032784212806261, 0.161440592119012, 0.0918604901195994, -0.172035353970668, -0.09236727461758, 0.227541309335466, -0.0271064288927961, -0.240067863523876, 0.328655132963512, -0.253548479406976, 0.138539534773853, -0.0156245814715821, -0.0270625689557795, -0.0347792443504004, -0.0404017567197903, -0.043697667166768, -0.0437111427196691, -0.0390673686861722, -0.028326560107256, -0.010565475144516, 0.013772667364788, 0.0418872839957572, 0.0678809552638599, 0.0829774305641306, 0.0775707382850504, 0.0456337753079536, -0.00926033705716051, -0.0693097821487021, -0.104484426238643, -0.0859914975851302, -0.00968602943656027, 0.0848094234024644, 0.124292030854912, 0.0568361406765031, -0.0775036208422535, -0.14328385339256, -0.0388909671650906, 0.133989460979655, 0.12457267377836, -0.0987442790957019, -0.167208797651234, 0.104358774971947, 0.174169756894794, -0.192988043926538, -0.0675499881391024, 0.292924218390552, -0.320229001699121, 0.220460802618907, -3.56876484224812e-247, -0.00102512171907325, -0.00397027934580623, -0.00936120581288402, -0.0174842548136789, -0.0282482737985146, -0.0409466739345152, -0.0539764105269099, -0.0646250288604565, -0.0691180519428952, -0.0631892187521977, -0.0434208536623169, -0.00938726162410613, 0.0338873112206848, 0.0744520416822103, 0.0949769840587381, 0.079225149881395, 0.0231917706155978, -0.0540102560187259, -0.109705057138822, -0.0980091333053639, -0.00902172836722212, 0.10019832246793, 0.128125147284499, 0.0227674513033914, -0.124017203072707, -0.127456847253707, 0.0545391846649376, 0.173668358008991, -0.000360570458808516, -0.198481806067928, 0.00720839003565937, 0.226568818044302, -0.124871489843427, -0.165127357662195, 0.334858788046061, -0.307587295788329, 0.0156133420594507, 0.0270431017229207, 0.034766373915688, 0.0404445609587174, 0.0439090301140181, 0.0442936954731884, 0.0403204276495454, 0.0306073443258258, 0.0141701605336169, -0.00883510287837906, -0.036214735482857, -0.0629865042940989, -0.0813448500391085, -0.0820969973836211, -0.0581516827368856, -0.00976059608569868, 0.0505321497155697, 0.0973336791519515, 0.10094561614332, 0.0461552550123792, -0.0468273060021331, -0.118463920594529, -0.102684306114828, 0.00923408976789581, 0.126927948710023, 0.117322877297659, -0.0401157840971917, -0.163498927546721, -0.0543390434493087, 0.161513908122053, 0.107877709923145, -0.174902751068176, -0.104781579279544, 0.239340770854956, -0.0199898409911556, -0.273771110171716, 0.378199707791041, 4.806803929172e-247, -0.00100301525933643, -0.0038846613954019, -0.00916198569518429, -0.0171282585169644, -0.0277315881064119, -0.0403601345779541, -0.0535821673434292, -0.0649335179590927, -0.0709141390251499, -0.0674204964867083, -0.0508440136545685, -0.0199170658738171, 0.0220082942861674, 0.0652219155871469, 0.0941356074446415, 0.0918229864295389, 0.0491724067029432, -0.0242028095673985, -0.0948764061998307, -0.116313450640292, -0.0593141927772445, 0.0518523589046972, 0.133009616761822, 0.0967117723101281, -0.0478413227213283, -0.1552193534347, -0.074460196512919, 0.124545225480018, 0.153632001496229, -0.0816740900925652, -0.194770271788664, 0.0876607868143915, 0.213209512208533, -0.190006382694524, -0.142777956445872, 0.41409049000431, 0.016158061323897, 0.0279865831648033, 0.0359918361784229, 0.0419297266439506, 0.0456915215481946, 0.0464737400277126, 0.0430560401212055, 0.0340771395865148, 0.0184767935622162, -0.00384150449640991, -0.0311828947240626, -0.0593003459727688, -0.0811635327666808, -0.0879118121200237, -0.0715910312175267, -0.029684010814238, 0.0298810740594332, 0.0863704869157036, 0.110811187786383, 0.0801653038331425, -0.00245819689237728, -0.0949101543490569, -0.129652398971799, -0.0622989064948419, 0.0709383604730409, 0.150256292366713, 0.0725361046759152, -0.104873712556326, -0.164503961161396, 0.0126074835227332, 0.194727858242832, 0.0540776885031338, -0.21647730988049, -0.064807101653981, 0.274962731116647, -0.0240986643481792, -0.40930210016891, 4.29940933555997e-247, 0.000999163663334846, 0.00386974422823145, 0.0091295129928173, 0.0170839915512596, 0.0277196469225884, 0.0405086898644248, 0.0541668257274653, 0.0664409148658943, 0.0740622489469013, 0.0730607193905577, 0.0596519878121899, 0.0318125992046044, -0.00863524251721547, -0.0541317486899773, -0.091018281981426, -0.102392297082215, -0.0752273984283337, -0.0104917720835538, 0.0688793284996559, 0.12067672939288, 0.103742277118375, 0.01188174187477, -0.1015013055614, -0.142560006375526, -0.0534672827871095, 0.103555248110814, 0.162170210863898, 0.023103405643826, -0.167713894807243, -0.130487773447006, 0.130714117652445, 0.196490268028899, -0.119171993161488, -0.245143472206481, 0.18265418283104, 0.370342296243705, -0.0169006531462098, -0.0292727899303342, -0.0376595642290761, -0.0439372856982485, -0.0480638569797865, -0.0493007917390541, -0.0464869953075899, -0.038284767072823, -0.0235648593726108, -0.00198998091195924, 0.0252329374335801, 0.0546019584090563, 0.0799342989825471, 0.0928935960353601, 0.0850073449432555, 0.0513988094491129, -0.00448294639717783, -0.0673268547078672, -0.111164694563174, -0.108642822640497, -0.048097926754037, 0.0492289959229715, 0.127199301307649, 0.122562460415146, 0.0182487412685992, -0.116098362340792, -0.154328862376577, -0.0313254626827251, 0.145579017779526, 0.159171850983335, -0.0502884779438395, -0.215189645203365, -0.0468729711518036, 0.240232464691157, 0.12057382407948, -0.29792462892547, -0.310232730942204, -2.81341724219163e-247, 0.000992333663421933, 0.00384329175231414, 0.0090699090106022, 0.0169894566286642, 0.0276280222757123, 0.0405464281950781, 0.0546183633137168, 0.0678228452455812, 0.0771614124229296, 0.0788630422147326, 0.0690651419433117, 0.045103667372627, 0.00734282581037968, -0.0389189129011332, -0.0824625842403483, -0.107506256858374, -0.0986650468944125, -0.0494092863142507, 0.0288988381384461, 0.103675809560517, 0.130965166142549, 0.0809946939849946, -0.0309854867436997, -0.134261903482456, -0.13941635913022, -0.0178605504382919, 0.13788581633435, 0.164030084417374, -0.00208648897755389, -0.190295901234224, -0.139921711639977, 0.133118552291814, 0.242630689916179, -0.046567033359139, -0.353911477066066, -0.242444852963543, -0.0177120047598799, -0.0306780921480137, -0.0394823301977585, -0.046134376212485, -0.0506683513142565, -0.0524225953362923, -0.0503106086796104, -0.0430373457256695, -0.0294260648285456, -0.00891763540252929, 0.017765882810155, 0.0479089956738412, 0.0763113054487843, 0.0954514437585438, 0.0968206205483316, 0.0737706747928635, 0.0256219674496607, -0.0383110867091523, -0.0973925293296832, -0.124228954734873, -0.0968284050087884, -0.0154403560151809, 0.0856477122928086, 0.145437798461342, 0.110900779347191, -0.0137423676014385, -0.143185571282607, -0.157587231219329, -0.0151651690988819, 0.164188256149857, 0.175737616330352, -0.0357650282864061, -0.236184818721952, -0.133743769708777, 0.199432804047876, 0.353231066019312, 0.177531103552202, -1.41478559155626e-247, -0.000979455695966268, -0.00379341559882535, -0.00895507516347714, -0.0167917606297128, -0.0273697625187614, -0.0403431809834604, -0.0547556685542432, -0.0688432965693369, -0.0799247088929789, -0.084509019048798, -0.0787810231148805, -0.0595996484478637, -0.0260186414932035, 0.0189348660413875, 0.0668988088800176, 0.104325670510366, 0.115360109504937, 0.0880378349426003, 0.0226545474034842, -0.0615234880052467, -0.127307962500369, -0.133731858655478, -0.0621064771614211, 0.0591226906465053, 0.154057403385876, 0.141642616377577, 0.00754296600359896, -0.151932722223746, -0.183023847999168, -0.0226972691919893, 0.189922831574814, 0.210424244548932, -0.0320017905949335, -0.298507529632765, -0.310962644130615, -0.122051600246348, 0.0186264886496393, 0.0322620247077803, 0.041536997616596, 0.048611964788949, 0.0536082392231164, 0.0559536770355448, 0.0546529316080015, 0.0484753115812207, 0.0362237965105906, 0.0171474590515715, -0.00849458596588017, -0.0387874031197759, -0.0696175943433108, -0.0945551502019377, -0.105581250859773, -0.0950401569561012, -0.058867819897613, -0.000443913519324001, 0.0665890121490754, 0.119057541598996, 0.130916760205112, 0.0860520626480827, -0.00707111104200686, -0.108437921425276, -0.15833512549733, -0.111599422122476, 0.0204309603706028, 0.154350727204375, 0.178227648826262, 0.0460238822520313, -0.15033750331705, -0.226939255658505, -0.0776820558676822, 0.187859271854652, 0.330750525105174, 0.246819459719516, 0.0788122157388073, -5.52547635008292e-248, 0.000960347108203595, 0.00371940835665098, 0.00878327916509135, 0.0164872776015251, 0.0269375908101224, 0.0398844620758831, 0.0545504700129918, 0.069449607117221, 0.0822609575550383, 0.0898559464275524, 0.0886067387092937, 0.0751009675927543, 0.0473194727998288, 0.00616309039114486, -0.0430704968452884, -0.0900057543942263, -0.120296514921443, -0.11937455591835, -0.0788396725570985, -0.00378809595197557, 0.0827468496109204, 0.143292738152927, 0.14072916334786, 0.0615976217989057, -0.0641193537183748, -0.165141346259337, -0.165254789357351, -0.0424021412971794, 0.131814922783182, 0.221940371591755, 0.132274492515749, -0.0930296454651355, -0.286042262074055, -0.306190044432458, -0.178598431610767, -0.0477625566072432, -0.0197004241776276, -0.0341221356063094, -0.0439495445021536, -0.0515193193196801, -0.0570533953677394, -0.060083406582858, -0.0597235787726066, -0.0548325810766063, -0.0442268317089168, -0.0270059258683004, -0.00300919831414675, 0.0266303852536889, 0.0589486634268553, 0.088845218421707, 0.109323079565022, 0.11261334417448, 0.092393102909526, 0.0469220883193455, -0.0177321431449971, -0.0859840670473247, -0.134604369077902, -0.139737353161091, -0.0885198073403454, 0.0088531046414966, 0.114361060897587, 0.172831250295598, 0.139973384019729, 0.0166742076204815, -0.134525100569405, -0.214244256423459, -0.151079330291256, 0.0348619809603696, 0.230430728712657, 0.311660224547996, 0.247297713993638, 0.118462883187103, 0.0271179190766295, 1.65957086095535e-248, 0.000933757015467194, 0.00361642537030879, 0.00854301263590735, 0.0160539936788618, 0.0262940502977177, 0.0391106979760226, 0.0539117920565591, 0.0695068023459449, 0.0839758529872224, 0.0946363402866892, 0.0982015747657276, 0.0912350748191281, 0.0709760841401324, 0.0365209263813405, -0.00982714848701344, -0.0614691127855564, -0.107356188036935, -0.133544786721546, -0.126978265147211, -0.0812074974586158, -0.00240867900182016, 0.0874176640530708, 0.153798341680204, 0.161889003045505, 0.0950914781036496, -0.0282723563350104, -0.152838604464256, -0.207810695329164, -0.147632785109006, 0.0115495488148331, 0.190709494002083, 0.294590852224643, 0.278267587728965, 0.177677912259221, 0.0721731658332491, 0.0143848011694515, -0.0210169267288778, -0.0364023849133687, -0.0469059872812353, -0.0550773172709367, -0.0612565098430048, -0.065096278090638, -0.0658423325499806, -0.0624729424427724, -0.0538605087981141, -0.039012734770843, -0.0174180266808541, 0.0105134641240943, 0.0429714602712858, 0.0763563633705684, 0.105181186273579, 0.122519609560525, 0.121235870007682, 0.0960784089763298, 0.0463837467174011, -0.0213624594031066, -0.0923792783709745, -0.145361157295409, -0.158186587867862, -0.117072034513371, -0.0262788211542727, 0.0868112900913831, 0.176782796240294, 0.197519815666095, 0.126526079267764, -0.0168705011359398, -0.174259563658306, -0.276722538772085, -0.285048222226271, -0.213223538688843, -0.114581931581652, -0.0403676950786324, -0.00710160082482026, -3.72645193912897e-249, 0.000897210829649788, 0.00347488260127057, 0.00821154856433363, 0.0154486743941261, 0.0253664217520631, 0.0379084339981621, 0.0526709294880072, 0.0687721563947051, 0.0847296078583233, 0.0983892894756125, 0.10696966594183, 0.107295524866528, 0.096292346638524, 0.0717747122464075, 0.0334781812640129, -0.0158523022687447, -0.0696752605807944, -0.117674277436726, -0.147124520623427, -0.145921170581934, -0.107022833667429, -0.03318769695994, 0.0600319833869589, 0.145465462810973, 0.191292039480531, 0.173042001653515, 0.0864368987271094, -0.0458012573250261, -0.179392085838621, -0.266857464490819, -0.280766084433221, -0.227806773359419, -0.142963965888313, -0.0665193320402209, -0.0206704761535859, -0.00324602964385963, 0.0227289756121569, 0.0393677405642496, 0.0507486802659843, 0.0596923732404993, 0.066682230565182, 0.0715146263611296, 0.073594308772447, 0.0720571324379131, 0.0658861802519876, 0.0540750184244937, 0.0358636585688228, 0.01105326548735, -0.0196210089769519, -0.0541178653934802, -0.0888246723101849, -0.118551899218352, -0.136964010744137, -0.137594892961124, -0.115480476825651, -0.0692193734665653, -0.00296424785056712, 0.072459433232492, 0.140247440409084, 0.180819037462153, 0.177442749056434, 0.122954082163302, 0.0249621618542502, -0.0935492585081039, -0.200197102415806, -0.264777983844301, -0.271804980627983, -0.227318695944955, -0.155233870716085, -0.0844865159894555, -0.034714129353058, -0.00964005937540164, -0.00136412404828164, -6.00198394879859e-250, -0.000846100735597875, -0.00327693405818442, -0.00774659178219464, -0.0145910032960095, -0.0240201221653966, -0.036069018944674, -0.0505203312212195, -0.0668075014102398, -0.0839134044096444, -0.100288124194646, -0.11382079276181, -0.121910655512166, -0.121687053581644, -0.110419918996787, -0.0861336570632762, -0.0483814989130324, 0.00094559267789857, 0.0569991016387647, 0.111847427357468, 0.155165156804264, 0.175947452026823, 0.165181967542066, 0.118981125132388, 0.0412123005808271, -0.0556101240677009, -0.152209353104041, -0.227397921023518, -0.264716750172857, -0.258357778460052, -0.21546923986797, -0.153164219491432, -0.0912033378101882, -0.0439892708421028, -0.0161951648357758, -0.00406156950031693, -0.00052205974261346, 0.0251770888537903, 0.0436079970814408, 0.0562395388868488, 0.0662683356445095, 0.0743620904893358, 0.0804938420010892, 0.0842650659539393, 0.0850198711504268, 0.0819251473386442, 0.0740671323083254, 0.0605874779347755, 0.0408700555854596, 0.0147779130460136, -0.0170766139888933, -0.053067629654633, -0.0903618006958124, -0.124882731242902, -0.151532294037732, -0.164752277681476, -0.159454295955823, -0.132243208698788, -0.082718403930947, -0.0144929042931635, 0.0645145810179955, 0.142905939912045, 0.207886498155595, 0.248298052788564, 0.25778303322468, 0.236920710361538, 0.193283993034159, 0.139058542301576, 0.0869625979117361, 0.0461628470835525, 0.0200363767459431, 0.00668441146465915, 0.00152707735148577, 0.000179529578809016, 5.78102797040991e-251, -0.000770462341379683, -0.00298398781704348, -0.00705671687153899, -0.0133075962852389, -0.0219654321121463, -0.0331454421752771, -0.0468056952802872, -0.0626885809324277, -0.0802533956761402, -0.0986095300658181, -0.116464522655822, -0.132106500537822, -0.143444841756304, -0.14813414702909, -0.1438020212281, -0.128387793585778, -0.100575067866443, -0.0602662921897568, -0.00900721006891619, 0.0497660697023453, 0.110797830341506, 0.167513662974277, 0.212980053118781, 0.241213386751574, 0.248566402165309, 0.234792108927774, 0.203369231418494, 0.160830048386923, 0.115157678431017, 0.0737125709237796, 0.0414187404513472, 0.0199011682416405, 0.00786299344241518, 0.00239851324568416, 0.000502780016434526, 5.4388784826608e-05, 0.0293863029891708, 0.050898569823857, 0.0656716691504383, 0.077523435539555, 0.0873924734476001, 0.0954888349614663, 0.101672181162737, 0.105572032328471, 0.106646852921664, 0.104232564676287, 0.0976000743573229, 0.086031713688248, 0.0689213211113175, 0.0458976224326222, 0.0169641196820339, -0.0173593969242589, -0.0559176342415605, -0.0968549728878163, -0.137641494005464, -0.175223079643026, -0.206314249256411, -0.227824471456272, -0.237369216666731, -0.233773651271215, -0.217443698403115, -0.190474318720436, -0.156404059396952, -0.119611970311026, -0.0844705180299614, -0.0544752923136326, -0.0316158664033285, -0.0161922123588262, -0.00712156872273205, -0.0025844688573947, -0.000726171140603572, -0.000140538681202497, -1.40628632589601e-05, 7.29361388997722e-252, 0.00064015966152098, 0.00247932770798453, 0.00586552499410413, 0.0110749114509656, 0.0183298128838627, 0.0277976505736698, 0.039579379272676, 0.0536904829789369, 0.0700365392786359, 0.0883856908881803, 0.108341313332773, 0.129319447339869, 0.150536911355232, 0.171017067713962, 0.189620498923119, 0.205106767731559, 0.216230395528875, 0.221868803982039, 0.221172310206193, 0.213717228932812, 0.199634630226543, 0.179682241141056, 0.155228769497239, 0.128131357853053, 0.100508520246875, 0.074439534220274, 0.05164911478368, 0.033252323628151, 0.0196287181657856, 0.0104626458515926, 0.0049354721996354, 0.00200418818949578, 0.000672930551404682, 0.000175240278305318, 3.14795286267501e-05, 2.92738550242375e-06, 0.0419929883820597, 0.0727339894393769, 0.0938884679074128, 0.111039650583824, 0.125763510176817, 0.138719592302122, 0.150202048321548, 0.160312440850567, 0.169031065586497, 0.176253553945383, 0.181814826160686, 0.185509452572003, 0.187112747384973, 0.186404721791803, 0.183197601758234, 0.177366347940186, 0.168880295316891, 0.157832669034297, 0.144463503122992, 0.129170699205986, 0.112504013687075, 0.0951380539549406, 0.0778231457497542, 0.0613171132755088, 0.0463059846663889, 0.033326213998073, 0.0227035932707985, 0.0145230410934974, 0.00863805719315237, 0.00471943111724355, 0.00233218609945907, 0.00102145011384433, 0.000385624348940016, 0.000120513842597142, 2.92402431232896e-05, 4.89818089824432e-06, 4.25076847121687e-07).finished();
}
} // namespace mrcpp
} // namespace detail
