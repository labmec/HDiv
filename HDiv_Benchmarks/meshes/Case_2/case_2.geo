h = 10*0.02835;
Geometry.Tolerance = 1e-05;
Mesh.Algorithm = 8;

// Define points
p0 = newp; Point(p0) = {0.5, 0.0, 1.0, h};
p1 = newp; Point(p1) = {0.5, 0.0, 0.0, h};
p2 = newp; Point(p2) = {0.5, 1.0, 0.0, h};
p3 = newp; Point(p3) = {0.5, 1.0, 1.0, h};
p4 = newp; Point(p4) = {0.0, 0.5, 0.0, h};
p5 = newp; Point(p5) = {1.0, 0.5, 0.0, h};
p6 = newp; Point(p6) = {1.0, 0.5, 1.0, h};
p7 = newp; Point(p7) = {0.0, 0.5, 1.0, h};
p8 = newp; Point(p8) = {0.0, 0.0, 0.5, h};
p9 = newp; Point(p9) = {1.0, 0.0, 0.5, h};
p10 = newp; Point(p10) = {1.0, 1.0, 0.5, h};
p11 = newp; Point(p11) = {0.0, 1.0, 0.5, h};
p12 = newp; Point(p12) = {0.75, 0.5, 1.0, h};
p13 = newp; Point(p13) = {0.75, 0.5, 0.5, h};
p14 = newp; Point(p14) = {0.75, 1.0, 0.5, h};
p15 = newp; Point(p15) = {0.75, 1.0, 1.0, h};
p16 = newp; Point(p16) = {0.5, 0.5, 0.75, h};
p17 = newp; Point(p17) = {1.0, 0.5, 0.75, h};
p18 = newp; Point(p18) = {1.0, 1.0, 0.75, h};
p19 = newp; Point(p19) = {0.5, 1.0, 0.75, h};
p20 = newp; Point(p20) = {0.5, 0.75, 0.5, h};
p21 = newp; Point(p21) = {1.0, 0.75, 0.5, h};
p22 = newp; Point(p22) = {1.0, 0.75, 1.0, h};
p23 = newp; Point(p23) = {0.5, 0.75, 1.0, h};
p24 = newp; Point(p24) = {0.5, 0.625, 0.5, h};
p25 = newp; Point(p25) = {0.75, 0.625, 0.5, h};
p26 = newp; Point(p26) = {0.75, 0.625, 0.75, h};
p27 = newp; Point(p27) = {0.5, 0.625, 0.75, h};
p28 = newp; Point(p28) = {0.625, 0.5, 0.75, h};
p29 = newp; Point(p29) = {0.625, 0.5, 0.5, h};
p30 = newp; Point(p30) = {0.625, 0.75, 0.5, h};
p31 = newp; Point(p31) = {0.625, 0.75, 0.75, h};
p32 = newp; Point(p32) = {0.5, 0.5, 0.625, h};
p33 = newp; Point(p33) = {0.75, 0.5, 0.625, h};
p34 = newp; Point(p34) = {0.75, 0.75, 0.625, h};
p35 = newp; Point(p35) = {0.5, 0.75, 0.625, h};
p36 = newp; Point(p36) = {0.0, 0.0, 1.0, h};
p37 = newp; Point(p37) = {0.0, 0.0, 0.0, h};
p38 = newp; Point(p38) = {0.0, 1.0, 0.0, h};
p39 = newp; Point(p39) = {0.0, 1.0, 1.0, h};
p40 = newp; Point(p40) = {1.0, 0.0, 1.0, h};
p41 = newp; Point(p41) = {1.0, 0.0, 0.0, h};
p42 = newp; Point(p42) = {1.0, 1.0, 0.0, h};
p43 = newp; Point(p43) = {1.0, 1.0, 1.0, h};
p44 = newp; Point(p44) = {0.5, 0.5, 0.0, h};
p45 = newp; Point(p45) = {0.5, 0.5, 1.0, h};
p46 = newp; Point(p46) = {0.5, 0.0, 0.5, h};
p47 = newp; Point(p47) = {0.5, 1.0, 0.5, h};
p48 = newp; Point(p48) = {1.0, 0.5, 0.5, h};
p49 = newp; Point(p49) = {0.0, 0.5, 0.5, h};
p50 = newp; Point(p50) = {0.75, 0.5, 0.75, h};
p51 = newp; Point(p51) = {0.75, 1.0, 0.75, h};
p52 = newp; Point(p52) = {0.75, 0.75, 0.5, h};
p53 = newp; Point(p53) = {0.75, 0.75, 1.0, h};
p54 = newp; Point(p54) = {1.0, 0.75, 0.75, h};
p55 = newp; Point(p55) = {0.5, 0.75, 0.75, h};
p56 = newp; Point(p56) = {0.625, 0.625, 0.5, h};
p57 = newp; Point(p57) = {0.625, 0.625, 0.75, h};
p58 = newp; Point(p58) = {0.75, 0.625, 0.625, h};
p59 = newp; Point(p59) = {0.5, 0.625, 0.625, h};
p60 = newp; Point(p60) = {0.625, 0.5, 0.625, h};
p61 = newp; Point(p61) = {0.625, 0.75, 0.625, h};
p62 = newp; Point(p62) = {0.5, 0.5, 0.5, h};
p63 = newp; Point(p63) = {0.75, 0.75, 0.75, h};
p64 = newp; Point(p64) = {0.625, 0.625, 0.625, h};
// End of point specification

pin100 = newp; Point(pin100) = {0.25, 0., 0., h};
pin010 = newp; Point(pin010) = {0., 0.25, 0., h};
pin001 = newp; Point(pin001) = {0., 0., 0.25, h};
pin110 = newp; Point(pin110) = {0.25, 0.25, 0., h};
pin101 = newp; Point(pin101) = {0.25, 0., 0.25, h};
pin011 = newp; Point(pin011) = {0., 0.25, 0.25, h};

pout100 = newp; Point(pout100) = {0.875, 1., 1., h};
pout010 = newp; Point(pout010) = {1., 0.875, 1., h};
pout001 = newp; Point(pout001) = {1., 1., 0.875, h};
pout110 = newp; Point(pout110) = {0.875, 0.875, 1., h};
pout101 = newp; Point(pout101) = {0.875, 1., 0.875, h};
pout011 = newp; Point(pout011) = {1., 0.875, 0.875, h};


// Define lines
frac_line_0= newl; Line(frac_line_0) = {p0, p36};
frac_line_1= newl; Line(frac_line_1) = {p0, p40};
frac_line_2= newl; Line(frac_line_2) = {p0, p45};
frac_line_3= newl; Line(frac_line_3) = {p0, p46};
frac_line_4= newl; Line(frac_line_4) = {p1, pin100};
frac_line_4ex= newl; Line(frac_line_4ex) = {pin100, p37};
frac_line_5= newl; Line(frac_line_5) = {p1, p41};
frac_line_6= newl; Line(frac_line_6) = {p1, p44};
frac_line_7= newl; Line(frac_line_7) = {p1, p46};
frac_line_8= newl; Line(frac_line_8) = {p2, p38};
frac_line_9= newl; Line(frac_line_9) = {p2, p42};
frac_line_10= newl; Line(frac_line_10) = {p2, p44};
frac_line_11= newl; Line(frac_line_11) = {p2, p47};
frac_line_12= newl; Line(frac_line_12) = {p3, p15};
frac_line_13= newl; Line(frac_line_13) = {p3, p19};
frac_line_14= newl; Line(frac_line_14) = {p3, p23};
frac_line_15= newl; Line(frac_line_15) = {p3, p39};
frac_line_16= newl; Line(frac_line_16) = {p4, pin010};
frac_line_16ex= newl; Line(frac_line_16ex) = {pin010, p37};
frac_line_17= newl; Line(frac_line_17) = {p4, p38};
frac_line_18= newl; Line(frac_line_18) = {p4, p44};
frac_line_19= newl; Line(frac_line_19) = {p4, p49};
frac_line_20= newl; Line(frac_line_20) = {p5, p41};
frac_line_21= newl; Line(frac_line_21) = {p5, p42};
frac_line_22= newl; Line(frac_line_22) = {p5, p44};
frac_line_23= newl; Line(frac_line_23) = {p5, p48};
frac_line_24= newl; Line(frac_line_24) = {p6, p12};
frac_line_25= newl; Line(frac_line_25) = {p6, p17};
frac_line_26= newl; Line(frac_line_26) = {p6, p22};
frac_line_27= newl; Line(frac_line_27) = {p6, p40};
frac_line_28= newl; Line(frac_line_28) = {p7, p36};
frac_line_29= newl; Line(frac_line_29) = {p7, p39};
frac_line_30= newl; Line(frac_line_30) = {p7, p45};
frac_line_31= newl; Line(frac_line_31) = {p7, p49};
frac_line_32= newl; Line(frac_line_32) = {p8, p36};
frac_line_33= newl; Line(frac_line_33) = {p8, pin001};
frac_line_33ex= newl; Line(frac_line_33ex) = {pin001, p37};
frac_line_34= newl; Line(frac_line_34) = {p8, p46};
frac_line_35= newl; Line(frac_line_35) = {p8, p49};
frac_line_36= newl; Line(frac_line_36) = {p9, p40};
frac_line_37= newl; Line(frac_line_37) = {p9, p41};
frac_line_38= newl; Line(frac_line_38) = {p9, p46};
frac_line_39= newl; Line(frac_line_39) = {p9, p48};
frac_line_40= newl; Line(frac_line_40) = {p10, p14};
frac_line_41= newl; Line(frac_line_41) = {p10, p18};
frac_line_42= newl; Line(frac_line_42) = {p10, p21};
frac_line_43= newl; Line(frac_line_43) = {p10, p42};
frac_line_44= newl; Line(frac_line_44) = {p11, p38};
frac_line_45= newl; Line(frac_line_45) = {p11, p39};
frac_line_46= newl; Line(frac_line_46) = {p11, p47};
frac_line_47= newl; Line(frac_line_47) = {p11, p49};
frac_line_48= newl; Line(frac_line_48) = {p12, p45};
frac_line_49= newl; Line(frac_line_49) = {p12, p50};
frac_line_50= newl; Line(frac_line_50) = {p12, p53};
frac_line_51= newl; Line(frac_line_51) = {p13, p25};
frac_line_52= newl; Line(frac_line_52) = {p13, p29};
frac_line_53= newl; Line(frac_line_53) = {p13, p33};
frac_line_54= newl; Line(frac_line_54) = {p13, p48};
frac_line_55= newl; Line(frac_line_55) = {p14, p47};
frac_line_56= newl; Line(frac_line_56) = {p14, p51};
frac_line_57= newl; Line(frac_line_57) = {p14, p52};
frac_line_58= newl; Line(frac_line_58) = {p15, pout100};
frac_line_58ex= newl; Line(frac_line_58ex) = {pout100, p43};
frac_line_59= newl; Line(frac_line_59) = {p15, p51};
frac_line_60= newl; Line(frac_line_60) = {p15, p53};
frac_line_61= newl; Line(frac_line_61) = {p16, p27};
frac_line_62= newl; Line(frac_line_62) = {p16, p28};
frac_line_63= newl; Line(frac_line_63) = {p16, p32};
frac_line_64= newl; Line(frac_line_64) = {p16, p45};
frac_line_65= newl; Line(frac_line_65) = {p17, p48};
frac_line_66= newl; Line(frac_line_66) = {p17, p50};
frac_line_67= newl; Line(frac_line_67) = {p17, p54};
frac_line_68= newl; Line(frac_line_68) = {p18, pout001};
frac_line_68ex= newl; Line(frac_line_68ex) = {pout001, p43};
frac_line_69= newl; Line(frac_line_69) = {p18, p51};
frac_line_70= newl; Line(frac_line_70) = {p18, p54};
frac_line_71= newl; Line(frac_line_71) = {p19, p47};
frac_line_72= newl; Line(frac_line_72) = {p19, p51};
frac_line_73= newl; Line(frac_line_73) = {p19, p55};
frac_line_74= newl; Line(frac_line_74) = {p20, p24};
frac_line_75= newl; Line(frac_line_75) = {p20, p30};
frac_line_76= newl; Line(frac_line_76) = {p20, p35};
frac_line_77= newl; Line(frac_line_77) = {p20, p47};
frac_line_78= newl; Line(frac_line_78) = {p21, p48};
frac_line_79= newl; Line(frac_line_79) = {p21, p52};
frac_line_80= newl; Line(frac_line_80) = {p21, p54};
frac_line_81= newl; Line(frac_line_81) = {p22, pout010};
frac_line_81ex= newl; Line(frac_line_81ex) = {pout010, p43};
frac_line_82= newl; Line(frac_line_82) = {p22, p53};
frac_line_83= newl; Line(frac_line_83) = {p22, p54};
frac_line_84= newl; Line(frac_line_84) = {p23, p45};
frac_line_85= newl; Line(frac_line_85) = {p23, p53};
frac_line_86= newl; Line(frac_line_86) = {p23, p55};
frac_line_87= newl; Line(frac_line_87) = {p24, p56};
frac_line_88= newl; Line(frac_line_88) = {p24, p59};
frac_line_89= newl; Line(frac_line_89) = {p24, p62};
frac_line_90= newl; Line(frac_line_90) = {p25, p52};
frac_line_91= newl; Line(frac_line_91) = {p25, p56};
frac_line_92= newl; Line(frac_line_92) = {p25, p58};
frac_line_93= newl; Line(frac_line_93) = {p26, p50};
frac_line_94= newl; Line(frac_line_94) = {p26, p57};
frac_line_95= newl; Line(frac_line_95) = {p26, p58};
frac_line_96= newl; Line(frac_line_96) = {p26, p63};
frac_line_97= newl; Line(frac_line_97) = {p27, p55};
frac_line_98= newl; Line(frac_line_98) = {p27, p57};
frac_line_99= newl; Line(frac_line_99) = {p27, p59};
frac_line_100= newl; Line(frac_line_100) = {p28, p50};
frac_line_101= newl; Line(frac_line_101) = {p28, p57};
frac_line_102= newl; Line(frac_line_102) = {p28, p60};
frac_line_103= newl; Line(frac_line_103) = {p29, p56};
frac_line_104= newl; Line(frac_line_104) = {p29, p60};
frac_line_105= newl; Line(frac_line_105) = {p29, p62};
frac_line_106= newl; Line(frac_line_106) = {p30, p52};
frac_line_107= newl; Line(frac_line_107) = {p30, p56};
frac_line_108= newl; Line(frac_line_108) = {p30, p61};
frac_line_109= newl; Line(frac_line_109) = {p31, p55};
frac_line_110= newl; Line(frac_line_110) = {p31, p57};
frac_line_111= newl; Line(frac_line_111) = {p31, p61};
frac_line_112= newl; Line(frac_line_112) = {p31, p63};
frac_line_113= newl; Line(frac_line_113) = {p32, p59};
frac_line_114= newl; Line(frac_line_114) = {p32, p60};
frac_line_115= newl; Line(frac_line_115) = {p32, p62};
frac_line_116= newl; Line(frac_line_116) = {p33, p50};
frac_line_117= newl; Line(frac_line_117) = {p33, p58};
frac_line_118= newl; Line(frac_line_118) = {p33, p60};
frac_line_119= newl; Line(frac_line_119) = {p34, p52};
frac_line_120= newl; Line(frac_line_120) = {p34, p58};
frac_line_121= newl; Line(frac_line_121) = {p34, p61};
frac_line_122= newl; Line(frac_line_122) = {p34, p63};
frac_line_123= newl; Line(frac_line_123) = {p35, p55};
frac_line_124= newl; Line(frac_line_124) = {p35, p59};
frac_line_125= newl; Line(frac_line_125) = {p35, p61};
frac_line_126= newl; Line(frac_line_126) = {p44, p62};
frac_line_127= newl; Line(frac_line_127) = {p46, p62};
frac_line_128= newl; Line(frac_line_128) = {p49, p62};
frac_line_129= newl; Line(frac_line_129) = {p51, p63};
frac_line_130= newl; Line(frac_line_130) = {p53, p63};
frac_line_131= newl; Line(frac_line_131) = {p54, p63};
frac_line_132= newl; Line(frac_line_132) = {p56, p64};
frac_line_133= newl; Line(frac_line_133) = {p57, p64};
frac_line_134= newl; Line(frac_line_134) = {p58, p64};
frac_line_135= newl; Line(frac_line_135) = {p59, p64};
frac_line_136= newl; Line(frac_line_136) = {p60, p64};
frac_line_137= newl; Line(frac_line_137) = {p61, p64};


in_line_xy1 = newl; Line(in_line_xy1) = {pin100, pin110};
in_line_xy2 = newl; Line(in_line_xy2) = {pin010, pin110};
in_line_xz1 = newl; Line(in_line_xz1) = {pin100, pin101};
in_line_xz2 = newl; Line(in_line_xz2) = {pin001, pin101};
in_line_yz1 = newl; Line(in_line_yz1) = {pin010, pin011};
in_line_yz2 = newl; Line(in_line_yz2) = {pin001, pin011};


out_line_xy1 = newl; Line(out_line_xy1) = {pout100, pout110};
out_line_xy2 = newl; Line(out_line_xy2) = {pout010, pout110};
out_line_xz1 = newl; Line(out_line_xz1) = {pout100, pout101};
out_line_xz2 = newl; Line(out_line_xz2) = {pout001, pout101};
out_line_yz1 = newl; Line(out_line_yz1) = {pout010, pout011};
out_line_yz2 = newl; Line(out_line_yz2) = {pout001, pout011};

// End of line specification

// Start domain specification
frac_loop_9 = newll;
Line Loop(frac_loop_9) = { frac_line_16, frac_line_16ex, -frac_line_33, -frac_line_33ex, frac_line_32, -frac_line_28, frac_line_29, -frac_line_45, frac_line_44, -frac_line_17};
auxiliary_9 = news; Plane Surface(auxiliary_9) = {frac_loop_9};

Line{frac_line_19} In Surface{auxiliary_9};
Line{frac_line_31} In Surface{auxiliary_9};
Line{frac_line_35} In Surface{auxiliary_9};
Line{frac_line_47} In Surface{auxiliary_9};

frac_loop_10 = newll;
Line Loop(frac_loop_10) = { frac_line_20, -frac_line_37, frac_line_36, -frac_line_27, frac_line_26, frac_line_81, frac_line_81ex, -frac_line_68, -frac_line_68ex, -frac_line_41, frac_line_43, -frac_line_21};
auxiliary_10 = news; Plane Surface(auxiliary_10) = {frac_loop_10};

Line{frac_line_23} In Surface{auxiliary_10};
Line{frac_line_25} In Surface{auxiliary_10};
Line{frac_line_39} In Surface{auxiliary_10};
Line{frac_line_42} In Surface{auxiliary_10};
Line{frac_line_65} In Surface{auxiliary_10};
Line{frac_line_67} In Surface{auxiliary_10};
Line{frac_line_70} In Surface{auxiliary_10};
Line{frac_line_78} In Surface{auxiliary_10};
Line{frac_line_80} In Surface{auxiliary_10};
Line{frac_line_83} In Surface{auxiliary_10};

frac_loop_11 = newll;
Line Loop(frac_loop_11) = { frac_line_0, -frac_line_32, frac_line_33, frac_line_33ex, -frac_line_4, -frac_line_4ex, frac_line_5, -frac_line_37, frac_line_36, -frac_line_1};
auxiliary_11 = news; Plane Surface(auxiliary_11) = {frac_loop_11};

Line{frac_line_3} In Surface{auxiliary_11};
Line{frac_line_7} In Surface{auxiliary_11};
Line{frac_line_34} In Surface{auxiliary_11};
Line{frac_line_38} In Surface{auxiliary_11};

frac_loop_12 = newll;
Line Loop(frac_loop_12) = { frac_line_8, -frac_line_44, frac_line_45, -frac_line_15, frac_line_12, frac_line_58, frac_line_58ex, -frac_line_68, -frac_line_68ex, -frac_line_41, frac_line_43, -frac_line_9};
auxiliary_12 = news; Plane Surface(auxiliary_12) = {frac_loop_12};

Line{frac_line_11} In Surface{auxiliary_12};
Line{frac_line_13} In Surface{auxiliary_12};
Line{frac_line_40} In Surface{auxiliary_12};
Line{frac_line_46} In Surface{auxiliary_12};
Line{frac_line_55} In Surface{auxiliary_12};
Line{frac_line_56} In Surface{auxiliary_12};
Line{frac_line_59} In Surface{auxiliary_12};
Line{frac_line_69} In Surface{auxiliary_12};
Line{frac_line_71} In Surface{auxiliary_12};
Line{frac_line_72} In Surface{auxiliary_12};

frac_loop_13 = newll;
Line Loop(frac_loop_13) = { frac_line_4, frac_line_4ex, -frac_line_16, -frac_line_16ex, frac_line_17, -frac_line_8, frac_line_9, -frac_line_21, frac_line_20, -frac_line_5};
auxiliary_13 = news; Plane Surface(auxiliary_13) = {frac_loop_13};

Line{frac_line_6} In Surface{auxiliary_13};
Line{frac_line_10} In Surface{auxiliary_13};
Line{frac_line_18} In Surface{auxiliary_13};
Line{frac_line_22} In Surface{auxiliary_13};

frac_loop_14 = newll;
Line Loop(frac_loop_14) = { frac_line_0, -frac_line_28, frac_line_29, -frac_line_15, frac_line_12, frac_line_58, frac_line_58ex, -frac_line_81, -frac_line_81ex, -frac_line_26, frac_line_27, -frac_line_1};
auxiliary_14 = news; Plane Surface(auxiliary_14) = {frac_loop_14};

Line{frac_line_2} In Surface{auxiliary_14};
Line{frac_line_14} In Surface{auxiliary_14};
Line{frac_line_24} In Surface{auxiliary_14};
Line{frac_line_30} In Surface{auxiliary_14};
Line{frac_line_48} In Surface{auxiliary_14};
Line{frac_line_50} In Surface{auxiliary_14};
Line{frac_line_60} In Surface{auxiliary_14};
Line{frac_line_82} In Surface{auxiliary_14};
Line{frac_line_84} In Surface{auxiliary_14};
Line{frac_line_85} In Surface{auxiliary_14};

domain_loop = newsl;
Surface Loop(domain_loop) = {auxiliary_9,auxiliary_10,auxiliary_11,auxiliary_12,auxiliary_13,auxiliary_14};
Volume(1) = {domain_loop};
Physical Volume("DOMAIN") = {1};
// End of domain specification

// Start fracture specification
frac_loop_0 = newll;
Line Loop(frac_loop_0) = { frac_line_2, -frac_line_84, -frac_line_14, frac_line_13, frac_line_71, -frac_line_11, frac_line_10, -frac_line_6, frac_line_7, -frac_line_3};
fracture_0 = news; Plane Surface(fracture_0) = {frac_loop_0};
//Physical Surface("FRACTURE_0") = {fracture_0};
Surface{fracture_0} In Volume{1};

Line{frac_line_61} In Surface{fracture_0};
Line{frac_line_63} In Surface{fracture_0};
Line{frac_line_64} In Surface{fracture_0};
Line{frac_line_73} In Surface{fracture_0};
Line{frac_line_74} In Surface{fracture_0};
Line{frac_line_76} In Surface{fracture_0};
Line{frac_line_77} In Surface{fracture_0};
Line{frac_line_86} In Surface{fracture_0};
Line{frac_line_88} In Surface{fracture_0};
Line{frac_line_89} In Surface{fracture_0};
Line{frac_line_97} In Surface{fracture_0};
Line{frac_line_99} In Surface{fracture_0};
Line{frac_line_113} In Surface{fracture_0};
Line{frac_line_115} In Surface{fracture_0};
Line{frac_line_123} In Surface{fracture_0};
Line{frac_line_124} In Surface{fracture_0};
Line{frac_line_126} In Surface{fracture_0};
Line{frac_line_127} In Surface{fracture_0};

frac_loop_1 = newll;
Line Loop(frac_loop_1) = { frac_line_18, -frac_line_22, frac_line_23, -frac_line_65, -frac_line_25, frac_line_24, frac_line_48, -frac_line_30, frac_line_31, -frac_line_19};
fracture_1 = news; Plane Surface(fracture_1) = {frac_loop_1};
//Physical Surface("FRACTURE_1") = {fracture_1};
Surface{fracture_1} In Volume{1};

Line{frac_line_49} In Surface{fracture_1};
Line{frac_line_52} In Surface{fracture_1};
Line{frac_line_53} In Surface{fracture_1};
Line{frac_line_54} In Surface{fracture_1};
Line{frac_line_62} In Surface{fracture_1};
Line{frac_line_63} In Surface{fracture_1};
Line{frac_line_64} In Surface{fracture_1};
Line{frac_line_66} In Surface{fracture_1};
Line{frac_line_100} In Surface{fracture_1};
Line{frac_line_102} In Surface{fracture_1};
Line{frac_line_104} In Surface{fracture_1};
Line{frac_line_105} In Surface{fracture_1};
Line{frac_line_114} In Surface{fracture_1};
Line{frac_line_115} In Surface{fracture_1};
Line{frac_line_116} In Surface{fracture_1};
Line{frac_line_118} In Surface{fracture_1};
Line{frac_line_126} In Surface{fracture_1};
Line{frac_line_128} In Surface{fracture_1};

frac_loop_2 = newll;
Line Loop(frac_loop_2) = { frac_line_34, -frac_line_38, frac_line_39, -frac_line_78, -frac_line_42, frac_line_40, frac_line_55, -frac_line_46, frac_line_47, -frac_line_35};
fracture_2 = news; Plane Surface(fracture_2) = {frac_loop_2};
//Physical Surface("FRACTURE_2") = {fracture_2};
Surface{fracture_2} In Volume{1};

Line{frac_line_51} In Surface{fracture_2};
Line{frac_line_52} In Surface{fracture_2};
Line{frac_line_54} In Surface{fracture_2};
Line{frac_line_57} In Surface{fracture_2};
Line{frac_line_74} In Surface{fracture_2};
Line{frac_line_75} In Surface{fracture_2};
Line{frac_line_77} In Surface{fracture_2};
Line{frac_line_79} In Surface{fracture_2};
Line{frac_line_87} In Surface{fracture_2};
Line{frac_line_89} In Surface{fracture_2};
Line{frac_line_90} In Surface{fracture_2};
Line{frac_line_91} In Surface{fracture_2};
Line{frac_line_103} In Surface{fracture_2};
Line{frac_line_105} In Surface{fracture_2};
Line{frac_line_106} In Surface{fracture_2};
Line{frac_line_107} In Surface{fracture_2};
Line{frac_line_127} In Surface{fracture_2};
Line{frac_line_128} In Surface{fracture_2};

frac_loop_3 = newll;
Line Loop(frac_loop_3) = { frac_line_49, -frac_line_116, -frac_line_53, frac_line_51, frac_line_90, -frac_line_57, frac_line_56, -frac_line_59, frac_line_60, -frac_line_50};
fracture_3 = news; Plane Surface(fracture_3) = {frac_loop_3};
//Physical Surface("FRACTURE_3") = {fracture_3};
Surface{fracture_3} In Volume{1};

Line{frac_line_92} In Surface{fracture_3};
Line{frac_line_93} In Surface{fracture_3};
Line{frac_line_95} In Surface{fracture_3};
Line{frac_line_96} In Surface{fracture_3};
Line{frac_line_117} In Surface{fracture_3};
Line{frac_line_119} In Surface{fracture_3};
Line{frac_line_120} In Surface{fracture_3};
Line{frac_line_122} In Surface{fracture_3};
Line{frac_line_129} In Surface{fracture_3};
Line{frac_line_130} In Surface{fracture_3};

frac_loop_4 = newll;
Line Loop(frac_loop_4) = { frac_line_61, frac_line_97, -frac_line_73, frac_line_72, -frac_line_69, frac_line_70, -frac_line_67, frac_line_66, -frac_line_100, -frac_line_62};
fracture_4 = news; Plane Surface(fracture_4) = {frac_loop_4};
//Physical Surface("FRACTURE_4") = {fracture_4};
Surface{fracture_4} In Volume{1};

Line{frac_line_93} In Surface{fracture_4};
Line{frac_line_94} In Surface{fracture_4};
Line{frac_line_96} In Surface{fracture_4};
Line{frac_line_98} In Surface{fracture_4};
Line{frac_line_101} In Surface{fracture_4};
Line{frac_line_109} In Surface{fracture_4};
Line{frac_line_110} In Surface{fracture_4};
Line{frac_line_112} In Surface{fracture_4};
Line{frac_line_129} In Surface{fracture_4};
Line{frac_line_131} In Surface{fracture_4};

frac_loop_5 = newll;
Line Loop(frac_loop_5) = { frac_line_75, frac_line_106, -frac_line_79, frac_line_80, -frac_line_83, frac_line_82, -frac_line_85, frac_line_86, -frac_line_123, -frac_line_76};
fracture_5 = news; Plane Surface(fracture_5) = {frac_loop_5};
//Physical Surface("FRACTURE_5") = {fracture_5};
Surface{fracture_5} In Volume{1};

Line{frac_line_108} In Surface{fracture_5};
Line{frac_line_109} In Surface{fracture_5};
Line{frac_line_111} In Surface{fracture_5};
Line{frac_line_112} In Surface{fracture_5};
Line{frac_line_119} In Surface{fracture_5};
Line{frac_line_121} In Surface{fracture_5};
Line{frac_line_122} In Surface{fracture_5};
Line{frac_line_125} In Surface{fracture_5};
Line{frac_line_130} In Surface{fracture_5};
Line{frac_line_131} In Surface{fracture_5};

frac_loop_6 = newll;
Line Loop(frac_loop_6) = { frac_line_87, -frac_line_91, frac_line_92, -frac_line_95, frac_line_94, -frac_line_98, frac_line_99, -frac_line_88};
fracture_6 = news; Plane Surface(fracture_6) = {frac_loop_6};
//Physical Surface("FRACTURE_6") = {fracture_6};
Surface{fracture_6} In Volume{1};

Line{frac_line_132} In Surface{fracture_6};
Line{frac_line_133} In Surface{fracture_6};
Line{frac_line_134} In Surface{fracture_6};
Line{frac_line_135} In Surface{fracture_6};

frac_loop_7 = newll;
Line Loop(frac_loop_7) = { frac_line_101, -frac_line_110, frac_line_111, -frac_line_108, frac_line_107, -frac_line_103, frac_line_104, -frac_line_102};
fracture_7 = news; Plane Surface(fracture_7) = {frac_loop_7};
//Physical Surface("FRACTURE_7") = {fracture_7};
Surface{fracture_7} In Volume{1};

Line{frac_line_132} In Surface{fracture_7};
Line{frac_line_133} In Surface{fracture_7};
Line{frac_line_136} In Surface{fracture_7};
Line{frac_line_137} In Surface{fracture_7};

frac_loop_8 = newll;
Line Loop(frac_loop_8) = { frac_line_113, -frac_line_124, frac_line_125, -frac_line_121, frac_line_120, -frac_line_117, frac_line_118, -frac_line_114};
fracture_8 = news; Plane Surface(fracture_8) = {frac_loop_8};
//Physical Surface("FRACTURE_8") = {fracture_8};
Surface{fracture_8} In Volume{1};

Line{frac_line_134} In Surface{fracture_8};
Line{frac_line_135} In Surface{fracture_8};
Line{frac_line_136} In Surface{fracture_8};
Line{frac_line_137} In Surface{fracture_8};

Line{in_line_xy1} In Surface{auxiliary_13};
Line{in_line_xy2} In Surface{auxiliary_13};
Line{in_line_xz1} In Surface{auxiliary_11};
Line{in_line_xz2} In Surface{auxiliary_11};
Line{in_line_yz1} In Surface{auxiliary_9};
Line{in_line_yz2} In Surface{auxiliary_9};

Line{out_line_xy1} In Surface{auxiliary_14};
Line{out_line_xy2} In Surface{auxiliary_14};
Line{out_line_xz1} In Surface{auxiliary_12};
Line{out_line_xz2} In Surface{auxiliary_12};
Line{out_line_yz1} In Surface{auxiliary_10};
Line{out_line_yz2} In Surface{auxiliary_10};
// End of fracture specification

Line Loop(187) = {19, -37, 150, -149};
//+
Plane Surface(188) = {187};
//+
Line Loop(188) = {146, -145, 6, -19};
//+
Plane Surface(189) = {188};
//+
Line Loop(189) = {147, -148, 37, -6};
//+
Plane Surface(190) = {189};
//+
Line Loop(190) = {63, -74, 154, -153};
//+
Plane Surface(191) = {190};
//+
Line Loop(191) = {152, -151, 63, -88};
//+
Plane Surface(192) = {191};
//+
Line Loop(192) = {155, -156, 74, -88};
//+
Plane Surface(193) = {192};
//+
Physical Surface("bc_inlet") = {190, 188, 189};
//+
Physical Surface("bc_outlet") = {179, 191, 193, 192};
//+
Line Loop(193) = {69, -3, 4, 134, -122, -68};
//+
Plane Surface(194) = {193};
//+
Line Loop(194) = {133, -134, -9, 8};
//+
Plane Surface(195) = {194};
//+
Line Loop(195) = {3, -52, -27, 30, -2};
//+
Plane Surface(196) = {195};
//+
Line Loop(196) = {30, -40, 43, -70, -28};
//+
Plane Surface(197) = {196};
//+
Line Loop(197) = {43, -26, 23, -41};
//+
Plane Surface(198) = {197};
//+
Line Loop(198) = {25, -8, 7, -23};
//+
Plane Surface(199) = {198};
//+
Line Loop(199) = {41, -7, 9, -42};
//+
Plane Surface(200) = {199};
//+
Line Loop(200) = {40, -2, 4, -42};
//+
Plane Surface(201) = {200};
//+
Line Loop(201) = {25, 133, -112, -56, 58, -26};
//+
Plane Surface(202) = {201};
//+
Line Loop(202) = {52, -69, 68, 122, -112, -56, 58, -70, -28, 27};
//+
Plane Surface(203) = {202};
//+
Line Loop(203) = {90, -72, -28, 29};
//+
Plane Surface(204) = {203};
//+
Line Loop(204) = {89, -54, -27, 29};
//+
Plane Surface(205) = {204};
//+
Line Loop(205) = {53, -100, 103, -137, -54};
//+
Plane Surface(206) = {205};
//+
Line Loop(206) = {138, -103, 100, -71, 72};
//+
Plane Surface(207) = {206};
//+
Line Loop(207) = {43, -58, 56, 112, -134, -42};
//+
Plane Surface(208) = {207};
//+
Line Loop(208) = {72, -86, 84, -70};
//+
Plane Surface(209) = {208};
//+
Line Loop(209) = {28, 71, -107, -67, 69, -52, -27};
//+
Line Loop(210) = {28, 71, -53, -27};
//+
Surface(210) = {210};
//+
Line Loop(211) = {53, -107, -67, 69, -52};
//+
Line Loop(212) = {84, -58, 55, 97, -85};
//+
Plane Surface(211) = {212};
//+
Line Loop(213) = {90, 138, -137, -89};
//+
Plane Surface(212) = {213};
//+
Line Loop(214) = {129, -138, -86, 85, -126};
//+
Plane Surface(213) = {214};
//+
Line Loop(215) = {129, -103, 100, -123, -57, 55, 97, -126};
//+
Plane Surface(214) = {215};
//+
Line Loop(216) = {71, -123, -57, 58, -70};
//+
Plane Surface(215) = {216};
//+
Line Loop(217) = {56, 110, -98, -55};
//+
Plane Surface(216) = {217};
//+
Line Loop(218) = {99, 141, -139, -98};
//+
Plane Surface(217) = {218};
//+
Line Loop(219) = {110, 139, -143, -111};
//+
Plane Surface(218) = {219};
//+
Line Loop(220) = {56, 111, -125, -57};
//+
Plane Surface(219) = {220};
//+
Line Loop(221) = {109, -125, 123, -107};
//+
Plane Surface(220) = {221};
//+
Line Loop(222) = {108, 140, -143, -109};
//+
Plane Surface(221) = {222};
//+
Line Loop(223) = {101, -108, 107, -100};
//+
Plane Surface(222) = {223};
//+
Line Loop(224) = {102, 141, -140, -101};
//+
Plane Surface(223) = {224};
//+
Line Loop(225) = {99, -124, -57, 55};
//+
Plane Surface(224) = {225};
//+
Line Loop(226) = {124, -102, 100, -123};
//+
Plane Surface(225) = {226};
//+
Line Loop(227) = {141, -143, -125, 124};
//+
Plane Surface(226) = {227};
//+
Plane Surface(227) = {202};
//+
Physical Surface("Region_1") = {196, 197, 201, 194, 208, 203, 198, 200, 199, 195, 202, 204, 209, 205, 212, 206, 214, 213, 211, 210, 215, 207, 217, 223, 218, 221, 222, 220, 225, 219, 216, 226, 224};
//+
Line Loop(228) = {8, -21, 18, 146, -145, -5};
//+
Plane Surface(228) = {228};
//+
Line Loop(229) = {12, -21, 20, -10};
//+
Plane Surface(229) = {229};
//+
Line Loop(230) = {18, 149, -150, -36, 39, -22};
//+
Plane Surface(230) = {230};
//+
Line Loop(231) = {20, -48, 51, -22};
//+
Plane Surface(231) = {231};
//+
Line Loop(232) = {39, -34, 31, -35};
//+
Plane Surface(232) = {232};
//+
Line Loop(233) = {34, -51, 49, -32};
//+
Plane Surface(233) = {233};
//+
Line Loop(234) = {50, -77, -15, 17, -49};
//+
Plane Surface(234) = {234};
//+
Line Loop(235) = {10, -48, 50, -13};
//+
Plane Surface(235) = {235};
//+
Line Loop(236) = {11, -47, 44, 59, -13};
//+
Plane Surface(236) = {236};
//+
Line Loop(237) = {77, -59, -44, 45, 73, 154, -153, -62, -14, 15};
//+
Plane Surface(237) = {237};
//+
Line Loop(238) = {45, 73, 156, -155, -87, -29, 28, 70, -84, -46};
//+
Plane Surface(238) = {238};
//+
Line Loop(239) = {24, -47, 46, 84, -26};
//+
Plane Surface(239) = {239};
//+
Line Loop(240) = {25, -12, 11, -24};
//+
Plane Surface(240) = {240};
//+
Line Loop(241) = {91, -33, 32, -17, 16};
//+
Plane Surface(241) = {241};
//+
Line Loop(242) = {3, -33, 31, -1};
//+
Plane Surface(242) = {242};
//+
Line Loop(243) = {1, -35, 38, -4};
//+
Plane Surface(243) = {243};
//+
Line Loop(244) = {38, -9, 5, 147, -148, -36};
//+
Plane Surface(244) = {244};
//+
Plane Surface(245) = {237};
//+
Line Loop(245) = {151, -152, -87, -29, 27, 52, -91, -16, 14, 62};
//+
Plane Surface(246) = {245};
//+
Physical Surface("BoundaryConditions_1") = {235, 236, 237, 234};
//+
Physical Surface("BoundaryConditions_2") = {246, 241, 242, 196};
//+
Physical Surface("BoundaryConditions_3") = {243, 244, 201, 200};
//+
Physical Surface("BoundaryConditions_4") = {228, 229, 199, 240};
//+
Physical Surface("BoundaryConditions_5") = {230, 232, 233, 231};
//+
Physical Surface("BoundaryConditions_6") = {239, 238, 197, 198};

//+

Physical Surface("FRACTURES") = {fracture_0,fracture_1,fracture_2,fracture_3,fracture_4,fracture_5,fracture_6,fracture_7,fracture_8};

//+
Line(157) = {48, 47};
//+
Line(158) = {46, 45};
//+
Line(159) = {49, 50};
//+
Physical Curve("FractureIntersections") = {159, 157, 158};
//+
Line(160) = {52, 51};
//+
Line(161) = {54, 53};
//+
Line(162) = {56, 55};
//+
Physical Curve("FractureIntersections") += {160, 162, 161};
//+
Line(163) = {62, 61};
//+
Line(164) = {58, 57};
//+
Line(165) = {60, 59};
//+
Physical Curve("FractureIntersections") += {163, 165, 164};
