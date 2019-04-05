// Gmsh project created on Sun Mar 31 19:13:00 2019
SetFactory("OpenCASCADE");

//Small fractures defintion 
p1 = newp; Point(p1) = {0.5, 0.75, 0.625};
p2 = newp; Point(p2) = {0.625, 0.75, 0.625};
p3 = newp; Point(p3) = {0.75, 0.75, 0.625};
p4 = newp; Point(p4) = {0.5, 0.625, 0.625};
p5 = newp; Point(p5) = {0.75, 0.625, 0.625};
p6 = newp; Point(p6) = {0.5, 0.5, 0.625};
p7 = newp; Point(p7) = {0.625, 0.5, 0.625};
p8 = newp; Point(p8) = {0.75, 0.5, 0.625};

p9 = newp; Point(p9) = {0.625, 0.75, 0.75};
p10 = newp; Point(p10) = {0.625, 0.75, 0.5};
p11 = newp; Point(p11) = {0.625, 0.625, 0.75};
p12 = newp; Point(p12) = {0.625, 0.625, 0.5};
p13 = newp; Point(p13) = {0.625, 0.5, 0.75};
p14 = newp; Point(p14) = {0.625, 0.5, 0.5};

p15 = newp; Point(p15) = {0.5, 0.625, 0.75};
p16 = newp; Point(p16) = {0.5, 0.625, 0.5};
p17 = newp; Point(p17) = {0.625, 0.625, 0.75};
p18 = newp; Point(p18) = {0.625, 0.625, 0.5};
p19 = newp; Point(p19) = {0.75, 0.625, 0.75};
p20 = newp; Point(p20) = {0.75, 0.625, 0.5};

//Center point small fractures
p21 = newp; Point(p21) = {0.625, 0.625, 0.625};

l1 = newl; Line(l1) = {p1,p2};
l2 = newl; Line(l2) = {p2,p3};
l3 = newl; Line(l3) = {p4,p21};
l4 = newl; Line(l4) = {p21,p5};
l5 = newl; Line(l5) = {p6,p7};
l6 = newl; Line(l6) = {p7,p8};

l7 = newl; Line(l7) = {p9,p2};
l8 = newl; Line(l8) = {p2,p10};
l9 = newl; Line(l9) = {p11,p21};
l10 = newl; Line(l10) = {p21,p12};
l11 = newl; Line(l11) = {p7,p13};
l12 = newl; Line(l12) = {p7,p14};

l13 = newl; Line(l13) = {p4,p15};
l14 = newl; Line(l14) = {p4,p16};
l15 = newl; Line(l15) = {p5,p19};
l16 = newl; Line(l16) = {p5,p20};

l17 = newl; Line(l17) = {p1,p4};
l18 = newl; Line(l18) = {p4,p6};
l17 = newl; Line(l17) = {p2,p21};
l18 = newl; Line(l18) = {p21,p7};
l19 = newl; Line(l19) = {p3,p5};
l20 = newl; Line(l20) = {p5,p8};

l21 = newl; Line(l21) = {p9,p11};
l22 = newl; Line(l22) = {p11,p13};
l23 = newl; Line(l23) = {p10,p12};
l24 = newl; Line(l24) = {p12,p14};

l25 = newl; Line(l25) = {p11,p15};
l26 = newl; Line(l26) = {p11,p19};
l27 = newl; Line(l27) = {p12,p16};
l28 = newl; Line(l28) = {p12,p20};

//Small fractures surfaces 
frac_1_1_29 = newll; Line Loop(frac_1_1_29) = {1, 19, -3, -17};
frac_1_2_30 = newll; Line Loop(frac_1_2_30) = {2, 21, -4, -19};
frac_1_3_31 = newll; Line Loop(frac_1_3_31) = {3, 20, -5, -18};
frac_1_4_32 = newll; Line Loop(frac_1_4_32) = {20, 6, -22, -4};

su_frac_1_1 = news; Plane Surface(su_frac_1_1) = {frac_1_1_29};
su_frac_1_2 = news; Plane Surface(su_frac_1_2) =
{frac_1_2_30};
su_frac_1_3 = news; Plane Surface(su_frac_1_3) = {frac_1_3_31};
su_frac_1_4 = news; Plane Surface(su_frac_1_4) = {frac_1_4_32};

frac_2_1_33 = newll; Line Loop(frac_2_1_33) = {7, 19, -9, -23};
frac_2_2_34 = newll; Line Loop(frac_2_2_34) = {8, 25, -10, -19};
frac_2_3_35 = newll; Line Loop(frac_2_3_35) = {9, 20, 11, -24};
frac_2_4_36 = newll; Line Loop(frac_2_4_36) = {10, 26, -12, -20};

su_frac_2_1 = news; Plane Surface(su_frac_2_1) = {frac_2_1_33};
su_frac_2_2 = news; Plane Surface(su_frac_2_2) =
{frac_2_2_34};
su_frac_2_3 = news; Plane Surface(su_frac_2_3) = {frac_2_3_35};
su_frac_2_4 = news; Plane Surface(su_frac_2_4) = {frac_2_4_36};

frac_3_1_37 = newll; Line Loop(frac_3_1_37) = {13, -27, 9, -3};
frac_3_2_38 = newll; Line Loop(frac_3_2_38) = {14, -29, -10, -3};
frac_3_3_39 = newll; Line Loop(frac_3_3_39) = {9, 4, 15, -28};
frac_3_4_40 = newll; Line Loop(frac_3_4_40) = {10, 30, -16, -4};

su_frac_3_1 = news; Plane Surface(su_frac_3_1) = {frac_3_1_37};
su_frac_3_2 = news; Plane Surface(su_frac_3_2) =
{frac_3_2_38};
su_frac_3_3 = news; Plane Surface(su_frac_3_3) = {frac_3_3_39};
su_frac_3_4 = news; Plane Surface(su_frac_3_4) = {frac_3_4_40};


//Medium fractures definition
p22 = newp; Point(p22) = {1.0, 1.0, 0.75};
p23 = newp; Point(p23) = {0.75, 1.0, 0.75};
p24 = newp; Point(p24) = {0.5, 1.0, 0.75};
p25 = newp; Point(p25) = {1.0, 0.75, 0.75};
p26 = newp; Point(p26) = {0.5, 0.75, 0.75};
p27 = newp; Point(p27) = {1.0, 0.5, 0.75};
p28 = newp; Point(p28) = {0.75, 0.5, 0.75};
p29 = newp; Point(p29) = {0.5, 0.5, 0.75};

p30 = newp; Point(p30) = {0.75, 1.0, 0.5};
p31 = newp; Point(p31) = {0.75, 1.0, 1.0};
p32 = newp; Point(p32) = {0.75, 0.75, 0.5};
p33 = newp; Point(p33) = {0.75, 0.75, 1.0};
p34 = newp; Point(p34) = {0.75, 0.5, 0.5};
p35 = newp; Point(p35) = {0.75, 0.5, 1.0};

p36 = newp; Point(p36) = {1.0, 0.75, 0.5};
p37 = newp; Point(p37) = {1.0, 0.75, 1.0};
p38 = newp; Point(p38) = {0.5, 0.75, 0.5};
p39 = newp; Point(p39) = {0.5, 0.75, 1.0};

//Center point medium fractures
p40 = newp; Point(p40) = {0.75, 0.75, 0.75};

l29 = newl; Line(l29) = {p23,p24};
l30 = newl; Line(l30) = {p22,p23};
l31 = newl; Line(l31) = {p9,p26};
l32 = newl; Line(l32) = {p9,p40};
l33 = newl; Line(l33) = {p25,p40};
l34 = newl; Line(l34) = {p13,p29};
l35 = newl; Line(l35) = {p13,p28};
l36 = newl; Line(l36) = {p27, p28};

l37 = newl; Line(l37) = {p23,p31};
l38 = newl; Line(l38) = {p23,p30};
l39 = newl; Line(l39) = {p33,p40};
l40 = newl; Line(l40) = {p3, p40};
l41 = newl; Line(l41) = {p3,p32};
l42 = newl; Line(l42) = {p28,p35};
l43 = newl; Line(l43) = {p8,p28};
l44 = newl; Line(l44) = {p8,p34};

l45 = newl; Line(l45) = {p26,p39};
l46 = newl; Line(l46) = {p1,p26};
l47 = newl; Line(l47) = {p1, p38};
l48 = newl; Line(l48) = {p25,p37};
l49 = newl; Line(l49) = {p25,p36};

l50 = newl; Line(l50) = {p24,p26};
l51 = newl; Line(l51) = {p15,p26};
l52 = newl; Line(l52) = {p15,p29};
l53 = newl; Line(l53) = {p23,p40};
l54 = newl; Line(l54) = {p40, p19};
l55 = newl; Line(l55) = {p19, p28};
l56 = newl; Line(l56) = {p22,p25};
l57 = newl; Line(l57) = {p25,p27};

l58 = newl; Line(l58) = {p31,p33};
l59 = newl; Line(l59) = {p33,p35};
l60 = newl; Line(l60) = {p30,p32};
l61 = newl; Line(l61) = {p32,p20};
l62 = newl; Line(l62) = {p20,p34};

l63 = newl; Line(l63) = {p33,p39};
l64 = newl; Line(l64) = {p33,p37};
l65 = newl; Line(l65) = {p10,p38};
l63 = newl; Line(l63) = {p10,p32};
l64 = newl; Line(l64) = {p32,p36};


//Medium fractures surfaces 
frac_5_1_65 = newll; Line Loop(frac_5_1_65) = {55, 76, -57, 58, -79};
frac_5_2_66 = newll; Line Loop(frac_5_2_66) = {56, 79, -59, -82};
frac_5_3_67 = newll; Line Loop(frac_5_3_67) = {57, -77, 78, -60, 61, -81, -80, -58};
frac_5_4_68 = newll; Line Loop(frac_5_4_68) = {59, 80, 81, -62, -83};

su_frac_5_1 = news; Plane Surface(su_frac_5_1) = {frac_5_1_65};
su_frac_5_2 = news; Plane Surface(su_frac_5_2) =
{frac_5_2_66};
su_frac_5_3 = news; Plane Surface(su_frac_5_3) = {frac_5_3_67};
su_frac_5_4 = news; Plane Surface(su_frac_5_4) = {frac_5_4_68};

frac_6_1_69 = newll; Line Loop(frac_6_1_69) = {63, 84, 65, -79};
frac_6_2_70 = newll; Line Loop(frac_6_2_70) = {64, 86, -67, 66, -79};
frac_6_3_71 = newll; Line Loop(frac_6_3_71) = {65, 80, 81, 68, -85};
frac_6_4_72 = newll; Line Loop(frac_6_4_72) = {66, 80, 81, -69, 70, -88, -87, -67};

su_frac_6_1 = news; Plane Surface(su_frac_6_1) = {frac_6_1_69};
su_frac_6_2 = news; Plane Surface(su_frac_6_2) =
{frac_6_2_70};
su_frac_6_3 = news; Plane Surface(su_frac_6_3) = {frac_6_3_71};
su_frac_6_4 = news; Plane Surface(su_frac_6_4) = {frac_6_4_72};

frac_7_1_73 = newll; Line Loop(frac_7_1_73) = {71, -89, 65, -58, 57};
frac_7_2_74 = newll; Line Loop(frac_7_2_74) = {72, -57, 58, -66, 67, -92, 91, -73};
frac_7_3_75 = newll; Line Loop(frac_7_3_75) = {65, -59, 74, -90};
frac_7_4_76 = newll; Line Loop(frac_7_4_76) = {66, -59, 75, -93, -67};

su_frac_7_1 = news; Plane Surface(su_frac_7_1) = {frac_7_1_73};
su_frac_7_2 = news; Plane Surface(su_frac_7_2) =
{frac_7_2_74};
su_frac_7_3 = news; Plane Surface(su_frac_7_3) = {frac_7_3_75};
su_frac_7_4 = news; Plane Surface(su_frac_7_4) = {frac_7_4_76};


//Big fractures defintion
p41 = newp; Point(p41) = {1.0, 1.0, 0.5};
p42 = newp; Point(p42) = {0.5, 1.0, 0.5};
p43 = newp; Point(p43) = {0.0, 1.0, 0.5};
p44 = newp; Point(p44) = {1.0, 0.5, 0.5};
p45 = newp; Point(p45) = {0.0, 0.5, 0.5};
p46 = newp; Point(p46) = {1.0, 0.0, 0.5};
p47 = newp; Point(p47) = {0.5, 0.0, 0.5};
p48 = newp; Point(p48) = {0.0, 0.0, 0.5};

p49 = newp; Point(p49) = {0.5, 1.0, 0.0};
p50 = newp; Point(p50) = {0.5, 1.0, 1.0};
p51 = newp; Point(p51) = {0.5, 0.5, 0.0};
p52 = newp; Point(p52) = {0.5, 0.5, 1.0};
p53 = newp; Point(p53) = {0.5, 0.0, 0.0};
p54 = newp; Point(p54) = {0.5, 0.0, 1.0};

p55 = newp; Point(p55) = {1.0, 0.5, 0.0};
p56 = newp; Point(p56) = {1.0, 0.5, 1.0};
p57 = newp; Point(p57) = {0.0, 0.5, 0.0};
p58 = newp; Point(p58) = {0.0, 0.5, 1.0};

//Center point big fractures
p59 = newp; Point(p59) = {0.5, 0.5, 0.5};

l65 = newl; Line(l65) = {p30,p41};
l66 = newl; Line(l66) = {p30,p42};
l67 = newl; Line(l67) = {p42,p43};
l68 = newl; Line(l68) = {p34,p44};
l69 = newl; Line(l69) = {p14,p34};
l70 = newl; Line(l70) = {p14,p59};
l71 = newl; Line(l71) = {p45,p59};
l72 = newl; Line(l72) = {p46,p47};
l73 = newl; Line(l73) = {p47,p48};

l74 = newl; Line(l74) = {p42,p49};
l75 = newl; Line(l75) = {p24,p42};
l76 = newl; Line(l76) = {p24,p50};
l77 = newl; Line(l77) = {p51,p59};
l78 = newl; Line(l78) = {p6,p59};
l79 = newl; Line(l79) = {p6,p29};
l80 = newl; Line(l80) = {p29,p52};
l81 = newl; Line(l81) = {p47,p53};
l82 = newl; Line(l82) = {p47,p54};

l83 = newl; Line(l83) = {p44,p55};
l84 = newl; Line(l84) = {p27,p44};
l85 = newl; Line(l85) = {p27,p56};
l86 = newl; Line(l86) = {p45,p57};
l87 = newl; Line(l87) = {p45,p58};

l88 = newl; Line(l88) = {p36,p41};
l89 = newl; Line(l89) = {p36,p44};
l90 = newl; Line(l90) = {p44,p46};
l91 = newl; Line(l91) = {p38,p42};
l92 = newl; Line(l92) = {p16,p38};
l93 = newl; Line(l93) = {p16,p59};
l94 = newl; Line(l94) = {p47,p59};
l95 = newl; Line(l95) = {p43,p45};
l96 = newl; Line(l96) = {p45,p48};

l97 = newl; Line(l97) = {p49,p51};
l98 = newl; Line(l98) = {p51,p53};
l99 = newl; Line(l99) = {p39,p50};
l100 = newl; Line(l100) = {p39,p52};
l101 = newl; Line(l101) = {p52,p54};

l102 = newl; Line(l102) = {p51,p55};
l103 = newl; Line(l103) = {p51,p57};
l104 = newl; Line(l104) = {p35,p56};
l105 = newl; Line(l105) = {p35,p52};
l106 = newl; Line(l106) = {p52,p58};


frac_8_1_77 = newll; Line Loop(frac_8_1_77) = {118, -141, 142, -121, -122, 123, -146, 145, 144, -119};
frac_8_2_78 = newll; Line Loop(frac_8_2_78) = {120, 148, 124, -146, 145, 144};
frac_8_3_79 = newll; Line Loop(frac_8_3_79) = {121, 143, 125, 147, -123, 122};
frac_8_4_80 = newll; Line Loop(frac_8_4_80) = {124, -147, 126, -149};

su_frac_8_1 = news; Plane Surface(su_frac_8_1) = {frac_8_1_77};
su_frac_8_2 = news; Plane Surface(su_frac_8_2) =
{frac_8_2_78};
su_frac_8_3 = news; Plane Surface(su_frac_8_3) = {frac_8_3_79};
su_frac_8_4 = news; Plane Surface(su_frac_8_4) = {frac_8_4_80};

frac_9_1_81 = newll; Line Loop(frac_9_1_81) = {127, 150, 130, -146, 145, 144};
frac_9_2_82 = newll; Line Loop(frac_9_2_82) = {128, -144, -145, 146, -131, 132, 133, -153, 152, -129};
frac_9_3_83 = newll; Line Loop(frac_9_3_83) = {130, -147, 134, -151};
frac_9_4_84 = newll; Line Loop(frac_9_4_84) = {131, -147, 135, -154, -133, -132};

su_frac_9_1 = news; Plane Surface(su_frac_9_1) = {frac_9_1_81};
su_frac_9_2 = news; Plane Surface(su_frac_9_2) =
{frac_9_2_82};
su_frac_9_3 = news; Plane Surface(su_frac_9_3) = {frac_9_3_83};
su_frac_9_4 = news; Plane Surface(su_frac_9_4) = {frac_9_4_84};

frac_10_1_85 = newll; Line Loop(frac_10_1_85) = {136, -155, 130, -123, 122, 121};
frac_10_2_86 = newll; Line Loop(frac_10_2_86) = {137, -121, -122, 123, -131, 132, 133, -158, 157, -138};
frac_10_3_87 = newll; Line Loop(frac_10_3_87) = {130, -124, 139, -156};
frac_10_4_88 = newll; Line Loop(frac_10_4_88) = {131, -124, 140, -159, -133, -132};

su_frac_10_1 = news; Plane Surface(su_frac_10_1) = {frac_10_1_85};
su_frac_10_2 = news; Plane Surface(su_frac_10_2) =
{frac_10_2_86};
su_frac_10_3 = news; Plane Surface(su_frac_10_3) = {frac_10_3_87};
su_frac_10_4 = news; Plane Surface(su_frac_10_4) = {frac_10_4_88};





//Volume small block 1
Curve Loop(180) = {61, -81, -28, 24};
Plane Surface(184) = {180};
Curve Loop(181) = {69, -81, -15, 22};
Plane Surface(185) = {181};
Curve Loop(182) = {69, -61, -11, 6};
Plane Surface(186) = {182};
Surface Loop(1) = {185, 53, 184, 45, 186, 38};
Volume(1) = {1};

//Volume small block 2
Curve Loop(183) = {80, -28, -23, 58};
Plane Surface(187) = {183};
Curve Loop(184) = {66, -58, 7, 2};
Plane Surface(188) = {184};
Curve Loop(185) = {15, -80, -66, 21};
Plane Surface(189) = {185};
Surface Loop(2) = {189, 36, 188, 43, 187, 53};
Volume(2) = {2};

//Volume small block 3
Curve Loop(186) = {26, 123, -146, -29};
Plane Surface(190) = {186};
Curve Loop(187) = {60, -78, -27, 24};
Plane Surface(191) = {187};
Curve Loop(188) = {60, -132, 5, 11};
Plane Surface(192) = {188};
Curve Loop(189) = {78, -132, -18, 13};
Plane Surface(193) = {189};
Surface Loop(3) = {192, 37, 193, 51, 191, 45};
Volume(3) = {3};

//Volume small block 4
Curve Loop(190) = {23, 27, 77, -57};
Plane Surface(194) = {190};
Curve Loop(191) = {1, -7, 57, -72};
Plane Surface(195) = {191};
Curve Loop(192) = {17, 13, 77, -72};
Plane Surface(196) = {192};
Curve Loop(193) = {3, -19, -1, 17};
Plane Surface(197) = {193};
Surface Loop(4) = {196, 195, 194, 43, 51, 35};
Volume(4) = {4};

//Volume small block 5
Curve Loop(194) = {70, -88, -16, 22};
Plane Surface(198) = {194};
Curve Loop(195) = {70, -122, -12, 6};
Plane Surface(199) = {195};
Curve Loop(196) = {88, -122, -26, 30};
Plane Surface(200) = {196};
Surface Loop(5) = {199, 46, 198, 54, 200, 38};
Volume(5) = {5};

//Volume small block 6
Curve Loop(197) = {16, -87, -67, 21};
Plane Surface(201) = {197};
Curve Loop(198) = {92, -67, -2, 8};
Plane Surface(202) = {198};
Curve Loop(199) = {87, -30, -25, 92};
Plane Surface(203) = {199};
Surface Loop(6) = {201, 54, 203, 44, 202, 36};
Volume(6) = {6};

//Volume small block 7
Curve Loop(200) = {123, -131, 5, 12};
Plane Surface(204) = {200};
Curve Loop(201) = {131, -146, -14, 18};
Plane Surface(205) = {201};
Curve Loop(202) = {26, 123, -146, -29};
Plane Surface(206) = {202};
Surface Loop(7) = {190, 204, 46, 37, 52, 205};
Volume(7) = {7};

//Volume small block 8
Curve Loop(203) = {91, -73, 1, 8};
Plane Surface(207) = {203};
Curve Loop(204) = {17, 14, 145, -73};
Plane Surface(208) = {204};
Curve Loop(205) = {29, 145, -91, 25};
Plane Surface(209) = {205};
Surface Loop(8) = {209, 208, 207, 35, 44, 52};
Volume(8) = {8};





//Inlet definition
p60 = newp; Point(p60) = {0.25, 0., 0.};
p61 = newp; Point(p61) = {0., 0.25, 0.};
p62 = newp; Point(p62) = {0., 0., 0.25};
p63 = newp; Point(p63) = {0.25, 0.25, 0.};
p64 = newp; Point(p64) = {0.25, 0., 0.25};
p65 = newp; Point(p65) = {0., 0.25, 0.25};



//Outlet Definition
p66= newp; Point(p66) = {0.875, 1., 1.};
p67 = newp; Point(p67) = {1., 0.875, 1.};
p68 = newp; Point(p68) = {1., 1., 0.875};
p69 = newp; Point(p69) = {0.875, 0.875, 1.};
p70 = newp; Point(p70) = {0.875, 1., 0.875};
p71 = newp; Point(p71) = {1., 0.875, 0.875};



//Domain definition
p72 = newp; Point(p72) = {1.0, 0.0, 1.0};
p73 = newp; Point(p73) = {1.0, 1.0, 1.0};
p74 = newp; Point(p74) = {0.0, 0.0, 1.0};
p75 = newp; Point(p75) = {0.0, 1.0, 1.0};
p76 = newp; Point(p76) = {1.0, 0.0, 0.0};
p77 = newp; Point(p77) = {1.0, 1.0, 0.0};
p78 = newp; Point(p78) = {0.0, 0.0, 0.0};
p79 = newp; Point(p79) = {0.0, 1.0, 0.0};

//Inlet lines definition
l107 = newl; Line(l107) = {p60,p64};
l108 = newl; Line(l108) = {p60,p78};
l109 = newl; Line(l109) = {p62,p78};
l110 = newl; Line(l110) = {p62,p64};
l111 = newl; Line(l111) = {p60,p63};
l112 = newl; Line(l112) = {p61,p63};
l113 = newl; Line(l113) = {p61,p78};
l114 = newl; Line(l114) = {p61,p65};
l115 = newl; Line(l115) = {p62,p65};

inlet_iz = newll; Line Loop(inlet_iz) = {210, -213, 212, -211};
inlet_bo = newll; Line Loop(inlet_bo) = {216, -212, 218, -217};
inlet_fr = newll; Line Loop(inlet_fr) = {214, -215, 216, -211};

su_inlet_iz = news; Plane Surface(su_inlet_iz) = {inlet_iz};
su_inlet_bo = news; Plane Surface(su_inlet_bo) = {inlet_bo};
su_inlet_fr = news; Plane Surface(su_inlet_fr) = {inlet_fr};

//Outlet lines definition
l116 = newl; Line(l116) = {p67,p73};
l117 = newl; Line(l117) = {p66,p73};
l118 = newl; Line(l118) = {p66,p70};
l119 = newl; Line(l119) = {p70,p68};
l120 = newl; Line(l120) = {p68,p73};
l121 = newl; Line(l121) = {p68,p71};
l122 = newl; Line(l122) = {p67,p71};
l123 = newl; Line(l123) = {p67,p69};
l124 = newl; Line(l124) = {p66,p69};

outlet_ba = newll; Line Loop(outlet_ba) = {225, -226, 233, -232};
outlet_ri = newll; Line Loop(outlet_ri) = {229, -226, 227, 228};
outlet_to = newll; Line Loop(outlet_to) = {231, -230, 229, -225};

su_outlet_ba = news; Plane Surface(su_outlet_ba) = {outlet_ba};
su_outlet_ri = news; Plane Surface(su_outlet_ri) = {outlet_ri};
su_outlet_to = news; Plane Surface(su_outlet_to) = {outlet_to};

//Domain lines definition
l125 = newl; Line(l125) = {p72,p56};
l126 = newl; Line(l126) = {p56,p37};
l127 = newl; Line(l127) = {p37,p67};
l128 = newl; Line(l128) = {p74,p58};
l129 = newl; Line(l129) = {p58,p75};
l130 = newl; Line(l130) = {p57,p61};
l131 = newl; Line(l131) = {p57,p79};
l132 = newl; Line(l132) = {p55,p76};
l133 = newl; Line(l133) = {p55,p77};
l134 = newl; Line(l134) = {p46,p72};
l135 = newl; Line(l135) = {p46,p76};
l136 = newl; Line(l136) = {p48,p74};
l137 = newl; Line(l137) = {p48,p62};
l138 = newl; Line(l138) = {p68,p22};
l139 = newl; Line(l139) = {p22,p41};
l140 = newl; Line(l140) = {p41,p77};
l141 = newl; Line(l141) = {p75,p43};
l142 = newl; Line(l142) = {p43,p79};
l143 = newl; Line(l143) = {p54,p72};
l144 = newl; Line(l144) = {p54,p74};
l145 = newl; Line(l145) = {p53,p76};
l146 = newl; Line(l146) = {p53,p60};
l147 = newl; Line(l147) = {p31,p66};
l148 = newl; Line(l148) = {p31,p50};
l149 = newl; Line(l149) = {p50,p75};
l150 = newl; Line(l150) = {p49,p77};
l151 = newl; Line(l151) = {p49,p79};


//Medium block 1
Curve Loop(237) = {157, 241, -90, 85};
Plane Surface(240) = {237};
Curve Loop(238) = {138, 241, -74, 83};
Plane Surface(241) = {238};
Curve Loop(239) = {138, -157, -68, -62};
Plane Surface(242) = {239};
Surface Loop(9) = {241, 240, 116, 108, 242, 101};
Volume(9) = {9};

//Medium block 2
Curve Loop(240) = {242, 225, -226, -262, 84, 90};
Plane Surface(243) = {240};
Curve Loop(241) = {74, 242, 225, -229, 253, 82};
Plane Surface(244) = {241};
Curve Loop(242) = {253, 56, 63, 262, 226, -229};
Plane Surface(245) = {242};
Surface Loop(10) = {244, 243, 245, 99, 106, 116};
Volume(10) = {10};

//Medium block 3
Curve Loop(243) = {85, 158, -153, -89};
Plane Surface(246) = {243};
Curve Loop(244) = {153, -133, -78, 77, 71};
Plane Surface(247) = {244};
Curve Loop(245) = {68, 158, -133, -60, 61};
Plane Surface(248) = {245};
Surface Loop(11) = {108, 246, 114, 247, 248, 187, 194, 184, 191};
Volume(11) = {11};

//Medium block 4
Curve Loop(246) = {84, 89, 152, -263};
Plane Surface(249) = {246};
Curve Loop(247) = {63, 263, -129, -55};
Plane Surface(250) = {247};
Curve Loop(248) = {76, 71, 152, -129};
Plane Surface(251) = {248};
Surface Loop(12) = {251, 98, 250, 249, 106, 114};
Volume(12) = {12};

//Medium block 8
Curve Loop(249) = {72, -76, 128, -144, -73};
Plane Surface(252) = {249};
Curve Loop(250) = {64, 119, -128, -55};
Plane Surface(253) = {250};
Curve Loop(251) = {91, 144, -119, 86, -92};
Plane Surface(254) = {251};
Surface Loop(13) = {252, 254, 107, 253, 98, 207, 195, 188, 202};
Volume(13) = {13};

//Medium block 5
Curve Loop(252) = {137, -121, -70, 69, -62};
Plane Surface(255) = {252};
Curve Loop(253) = {142, -137, -83, 75};
Plane Surface(256) = {253};
Curve Loop(254) = {142, -121, -88, -87, 93};
Plane Surface(257) = {254};
Surface Loop(14) = {256, 117, 257, 255, 101, 198, 185, 189, 201};
Volume(14) = {14};

//Medium block 6
Curve Loop(255) = {75, 141, -254, 82};
Plane Surface(258) = {255};
Curve Loop(256) = {141, -118, 86, 93};
Plane Surface(259) = {256};
Curve Loop(257) = {254, -118, -64, -56};
Plane Surface(260) = {257};
Surface Loop(15) = {258, 260, 99, 107, 259, 117};
Volume(15) = {15};

//Big block 1
Curve Loop(258) = {240, -157, 158, 154, 258};
Plane Surface(261) = {258};
Curve Loop(259) = {249, -258, -135, -125};
Plane Surface(262) = {259};
Curve Loop(260) = {143, 249, 240, -138, 137};
Plane Surface(263) = {260};
Surface Loop(16) = {263, 166, 262, 261, 175, 242, 255, 199, 248, 186, 192, 204};
Volume(16) = {16};

//Big block 3
Curve Loop(261) = {154, 259, 243, -159};
Plane Surface(264) = {261};
Curve Loop(262) = {243, -140, 149, 251};
Plane Surface(265) = {262};
Curve Loop(263) = {126, 251, -259, -135};
Plane Surface(266) = {263};
Surface Loop(17) = {175, 264, 266, 167, 265, 183};
Volume(17) = {17};


//Big block 4
Curve Loop(264) = {153, 159, 244, -264, -152};
Plane Surface(267) = {264};
Curve Loop(265) = {264, 256, -120, -128, 129};
Plane Surface(268) = {265};
Curve Loop(266) = {244, 256, 148, 140};
Plane Surface(269) = {266};
Surface Loop(18) = {267, 269, 268, 165, 252, 251, 208, 196, 193, 205, 183, 247};
Volume(18) = {18};

//Big block 5
Curve Loop(267) = {250, -247, -136, 143};
Plane Surface(270) = {267};
Curve Loop(268) = {250, -260, -134, -125};
Plane Surface(271) = {268};
Curve Loop(269) = {151, 260, -247, -155};
Plane Surface(272) = {269};
Surface Loop(19) = {166, 270, 271, 272, 174, 180};
Volume(19) = {19};

//Big block 6 
Curve Loop(270) = {136, 248, -255, -141, 142};
Plane Surface(273) = {270};
Curve Loop(271) = {255, -265, -127, -119, 118};
Plane Surface(274) = {271};
Curve Loop(272) = {155, 248, -265, 150};
Plane Surface(275) = {272};
Surface Loop(20) = {275, 273, 274, 172, 180, 257, 259, 254, 200, 190, 203, 209};
Volume(20) = {20};

//Big block 7
Curve Loop(273) = {134, 261, 211, -212, -252, -126};
Plane Surface(276) = {273};
Curve Loop(274) = {252, 212, -216, -245, -139, 149};
Plane Surface(277) = {274};
Curve Loop(275) = {261, 211, -216, -245, -156, 151};
Plane Surface(278) = {275};
Surface Loop(21) = {276, 278, 277, 182, 167, 174};
Volume(21) = {21};

//Big block 8 
Curve Loop(276) = {127, 266, -257, -120};
Plane Surface(279) = {276};
Curve Loop(277) = {148, 139, 246, -257};
Plane Surface(280) = {277};
Curve Loop(278) = {156, 246, -266, 150};
Plane Surface(281) = {278};
Surface Loop(22) = {280, 281, 279, 172, 182, 165};
Volume(22) = {22};



//Domain minus inlet, outlet
Curve Loop(279) = {232, -233, -262, 84, 90, 242};
Plane Surface(282) = {279};
Curve Loop(280) = {56, 63, 262, 227, 228, 253};
Plane Surface(283) = {280};
Curve Loop(281) = {261, 214, -215, -245, -156, 151};
Plane Surface(284) = {281};

Curve Loop(282) = {82, 74, 242, 231, -230, 253};
Plane Surface(285) = {282};
Curve Loop(283) = {252, 213, -210, -261, -134, 126};
Plane Surface(286) = {283};
Curve Loop(284) = {149, 252, 218, -217, -245, -139};
Plane Surface(287) = {284};

//Physical especification
Physical Volume("Permeability_1") = {21, 22, 19, 20, 16, 10, 15, 14, 9, 13, 8, 7, 6, 5, 2, 1};

Physical Volume("Permeability_2") = {18, 17, 12, 11, 4, 3};


Physical Surface("Inlet_Surface") = {224, 222, 223};

Physical Surface("Outlet_Surface") = {237, 239, 238};

Physical Surface("Domain_Surface") = {272, 281, 275, 274, 279, 260, 253, 250, 268, 249, 267, 264, 261, 246, 240, 262, 271, 266, 265, 269, 280, 270, 273, 263, 256, 258, 241, 281, 282, 283, 284,285, 286, 287};

Physical Surface("Fractures") = {35, 36, 37, 38, 43, 44, 46, 45, 51, 52, 53, 54, 99, 98, 100, 101, 107, 106, 109, 108, 117, 116, 115, 114, 164, 165, 166, 167, 172, 173, 174, 175, 180, 181, 182, 183};

Physical Curve("Fracture_Intersection") = {20, 19, 3, 4, 59, 58, 79, 80, 123, 147, 146, 124, 145, 144, 122, 121, 57, 81, 130, 131, 132, 133, 10, 9, 67, 66, 65, 8, 1, 2, 7, 25, 26, 29, 30, 17, 18, 13, 14, 5, 6, 27, 28, 76, 77, 78, 60, 61, 62, 91, 92, 93, 73, 72, 71, 21, 22, 15, 16, 87, 88, 23, 24, 12, 11, 68, 69, 70};

//+
Physical Point("Points_Intersection") = {59, 40, 21};
