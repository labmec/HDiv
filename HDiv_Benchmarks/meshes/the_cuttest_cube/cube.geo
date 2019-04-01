// Gmsh project created on Sun Mar 31 19:13:00 2019
SetFactory("OpenCASCADE");

p1 = newp; Point(p1) = {0,0,0};
p2 = newp; Point(p2) = {0.5,0,0};
p3 = newp; Point(p3) = {1,0,0};
p4 = newp; Point(p4) = {0,0,0.5};
p5 = newp; Point(p5) = {0.5,0,0.5};
p6 = newp; Point(p6) = {1,0,0.5};
p7 = newp; Point(p7) = {0,0,1};
p8 = newp; Point(p8) = {0.5,0,1};
p9 = newp; Point(p9) = {1,0,1};

p10 = newp; Point(p10) = {0,0.5,0};
p11 = newp; Point(p11) = {0.5,0.5,0};
p12 = newp; Point(p12) = {1,0.5,0};
p13 = newp; Point(p13) = {0,0.5,0.5};
p14 = newp; Point(p14) = {0.5,0.5,0.5};
p15 = newp; Point(p15) = {1,0.5,0.5};
p16 = newp; Point(p16) = {0,0.5,1};
p17 = newp; Point(p17) = {0.5,0.5,1};
p18 = newp; Point(p18) = {1,0.5,1};

p19 = newp; Point(p19) = {0,1,0};
p20 = newp; Point(p20) = {0.5,1,0};
p21 = newp; Point(p21) = {1,1,0};
p22 = newp; Point(p22) = {0,1,0.5};
p23 = newp; Point(p23) = {0.5,1,0.5};
p24 = newp; Point(p24) = {1,1,0.5};
p25 = newp; Point(p25) = {0,1,1};
p26 = newp; Point(p26) = {0.5,1,1};
p27 = newp; Point(p27) = {1,1,1};

l1 = newl; Line(l1) = {p1,p2};
l2 = newl; Line(l2) = {p2,p3};
l3 = newl; Line(l3) = {p4,p5};
l4 = newl; Line(l4) = {p5,p6};
l5 = newl; Line(l5) = {p7,p8};
l6 = newl; Line(l6) = {p8,p9};

l7 = newl; Line(l7) = {p10,p11};
l8 = newl; Line(l8) = {p11,p12};
l9 = newl; Line(l9) = {p13,p14};
l10 = newl; Line(l10) = {p14,p15};
l11 = newl; Line(l11) = {p16,p17};
l12 = newl; Line(l12) = {p17,p18};

l13 = newl; Line(l13) = {p19,p20};
l14 = newl; Line(l14) = {p20,p21};
l15 = newl; Line(l15) = {p22,p23};
l16 = newl; Line(l16) = {p23,p24};
l17 = newl; Line(l17) = {p25,p26};
l18 = newl; Line(l18) = {p26,p27};

l19 = newl; Line(19) = {p1,p4};
l20 = newl; Line(20) = {p4,p7};
l21 = newl; Line(21) = {p2,p5};
l22 = newl; Line(22) = {p5,p8};
l23 = newl; Line(23) = {p3,p6};
l24 = newl; Line(24) = {p6,p9};

l25 = newl; Line(25) = {p10,p13};
l26 = newl; Line(26) = {p13,p16};
l27 = newl; Line(27) = {p11,p14};
l28 = newl; Line(28) = {p14,p17};
l29 = newl; Line(29) = {p12,p15};
l30 = newl; Line(30) = {p15,p18};

l31 = newl; Line(31) = {p19,p22};
l32 = newl; Line(32) = {p22,p25};
l33 = newl; Line(33) = {p20,p23};
l34 = newl; Line(34) = {p23,p26};
l35 = newl; Line(35) = {p21,p24};
l36 = newl; Line(36) = {p24,p27};

l37 = newl; Line(37) = {p1,p10};
l38 = newl; Line(38) = {p10,p19};
l39 = newl; Line(39) = {p2,p11};
l40 = newl; Line(40) = {p11,p20};
l41 = newl; Line(41) = {p3,p12};
l42 = newl; Line(42) = {p12,p21};

l43 = newl; Line(43) = {p4,p13};
l44 = newl; Line(44) = {p13,p22};
l45 = newl; Line(45) = {p5,p14};
l46 = newl; Line(46) = {p14,p23};
l47 = newl; Line(47) = {p6,p15};
l48 = newl; Line(48) = {p15,p24};

l49 = newl; Line(49) = {p7,p16};
l50 = newl; Line(50) = {p16,p25};
l51 = newl; Line(51) = {p8,p17};
l52 = newl; Line(52) = {p17,p26};
l53 = newl; Line(53) = {p9,p18};
l54 = newl; Line(54) = {p18,p27};



//+
Curve Loop(2) = {1, 39, -7, -37};
//+
Plane Surface(1) = {2};
//+
Curve Loop(3) = {2, 41, -8, -39};
//+
Plane Surface(2) = {3};
//+
Curve Loop(4) = {7, 40, -13, -38};
//+
Plane Surface(3) = {4};
//+
Curve Loop(5) = {8, 42, -14, -40};
//+
Plane Surface(4) = {5};
//+
Physical Surface("BC_Back") = {2, 1, 4, 3};
//+
Curve Loop(6) = {6, 53, -12, -51};
//+
Plane Surface(5) = {6};
//+
Curve Loop(7) = {5, 51, -11, -49};
//+
Plane Surface(6) = {7};
//+
Curve Loop(8) = {12, 54, -18, -52};
//+
Plane Surface(7) = {8};
//+
Curve Loop(9) = {11, 52, -17, -50};
//+
Plane Surface(8) = {9};
//+
Physical Surface("BC_Front") = {5, 6, 7, 8};
//+
Curve Loop(10) = {23, 47, -29, -41};
//+
Plane Surface(9) = {10};
//+
Curve Loop(11) = {24, 53, -30, -47};
//+
Plane Surface(10) = {11};
//+
Curve Loop(12) = {29, 48, -35, -42};
//+
Plane Surface(11) = {12};
//+
Curve Loop(13) = {30, 54, -36, -48};
//+
Plane Surface(12) = {13};
//+
Curve Loop(14) = {20, 49, -26, -43};
//+
Plane Surface(13) = {14};
//+
Curve Loop(15) = {19, 43, -25, -37};
//+
Plane Surface(14) = {15};
//+
Curve Loop(16) = {26, 50, -32, -44};
//+
Plane Surface(15) = {16};
//+
Curve Loop(17) = {25, 44, -31, -38};
//+
Plane Surface(16) = {17};
//+
Curve Loop(18) = {35, -16, -33, 14};
//+
Plane Surface(17) = {18};
//+
Curve Loop(19) = {36, -18, -34, 16};
//+
Plane Surface(18) = {19};
//+
Curve Loop(20) = {33, -15, -31, 13};
//+
Plane Surface(19) = {20};
//+
Curve Loop(21) = {34, -17, -32, 15};
//+
Plane Surface(20) = {21};
//+
Curve Loop(22) = {24, -6, -22, 4};
//+
Plane Surface(21) = {22};
//+
Curve Loop(23) = {23, -4, -21, 2};
//+
Plane Surface(22) = {23};
//+
Curve Loop(24) = {22, -5, -20, 3};
//+
Plane Surface(23) = {24};
//+
Curve Loop(25) = {21, -3, -19, 1};
//+
Plane Surface(24) = {25};
//+
Physical Surface("BC_Neumann") = {9, 10, 11, 12, 13, 14, 15, 16, 21, 23, 22, 24, 20, 18, 19, 17};
//+
Curve Loop(26) = {25, 9, -27, -7};
//+
Plane Surface(25) = {26};
//+
Curve Loop(27) = {8, 29, -10, -27};
//+
Plane Surface(26) = {27};
//+
Curve Loop(28) = {9, 28, -11, -26};
//+
Plane Surface(27) = {28};
//+
Curve Loop(29) = {10, 30, -12, -28};
//+
Plane Surface(28) = {29};
//+
Curve Loop(30) = {22, 51, -28, -45};
//+
Plane Surface(29) = {30};
//+
Curve Loop(31) = {21, 45, -27, -39};
//+
Plane Surface(30) = {31};
//+
Curve Loop(32) = {28, 52, -34, -46};
//+
Plane Surface(31) = {32};
//+
Curve Loop(33) = {27, 46, -33, -40};
//+
Plane Surface(32) = {33};
//+
Curve Loop(34) = {4, 47, -10, -45};
//+
Plane Surface(33) = {34};
//+
Curve Loop(35) = {3, 45, -9, -43};
//+
Plane Surface(34) = {35};
//+
Curve Loop(36) = {10, 48, -16, -46};
//+
Plane Surface(35) = {36};
//+
Curve Loop(37) = {9, 46, -15, -44};
//+
Plane Surface(36) = {37};



//+
Surface Loop(1) = {21, 10, 5, 29, 33, 28};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {22, 9, 2, 26, 33, 30};
//+
Volume(2) = {2};
//+
Surface Loop(3) = {23, 6, 13, 27, 34, 29};
//+
Volume(3) = {3};
//+
Surface Loop(4) = {24, 14, 1, 30, 25, 34};
//+
Volume(4) = {4};
//+
Surface Loop(5) = {7, 12, 18, 35, 31, 28};
//+
Volume(5) = {5};
//+
Surface Loop(6) = {17, 11, 4, 35, 32, 26};
//+
Volume(6) = {6};
//+
Surface Loop(7) = {19, 16, 3, 36, 25, 32};
//+
Volume(7) = {7};
//+
Surface Loop(8) = {20, 8, 15, 31, 27, 36};
//+
Volume(8) = {8};
//+
Physical Volume("RockMatrix") = {8, 3, 4, 7, 2, 6, 1, 5};

//+
Physical Surface("Fractures") = {33, 34, 35, 36, 26, 28, 25, 27, 32, 31, 30, 29};
Physical Curve("FracturesIntersections") = {45, 28, 9, 46, 27, 10};
Physical Point("CrossingIntresections") = {14};

Transfinite Line "*" = 2 Using Bump 1;
Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";//+

Coherence;
