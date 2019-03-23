// Gmsh project created on Fri Mar 22 17:07:30 2019
SetFactory("OpenCASCADE");

n=2;
Geometry.Tolerance = 1e-05;
Mesh.Algorithm = 6;

//Back side domain points
Point(1)={0,0,0};
Point(2)={1,0,0};
Point(3)={1,1,0};
Point(4)={0,1,0};

//Back side domain lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

//Back side domain surfaces
Curve Loop(1) = {3, 4, 1, 2};
Plane Surface(1) = {1};


//Right side domain points
Point(5)={1,0,1};
Point(6)={1,1,1};

//Right side domain lines
Line(5) = {2, 5};
Line(6) = {5, 6};
Line(7) = {6, 3};

//Right side domain surfaces
Curve Loop(2) = {2, -7, -6, -5};
Plane Surface(2) = {2};


//Front side domain points
Point(7)={0,0,1};
Point(8)={0,1,1};

//Front side domain lines
Line(8) = {6, 8};
Line(9) = {8, 7};
Line(10) = {7, 5};

//Front side domain surfaces
Curve Loop(3) = {6, 8, -9, -10};
Plane Surface(3) = {3};


//Left side domain points (NONE)

//Left side domain lines
Line(11) = {8, 4};
Line(12) = {1, 7};

//Left side domain surfaces
Curve Loop(4) = {-9, -11, -4, 12};
Plane Surface(4) = {4};


//Top side domain points (NONE)

//Top side domain lines (NONE)

//Top side domain surfaces
Curve Loop(6) = {3, -11, -8, 7};
Plane Surface(5) = {6};


//Bottom side domain points (NONE)

//Bottom side domain lines (NONE)

//Bottom side domain surfaces
Curve Loop(7) = {1, 5, -10, -12};
Plane Surface(6) = {7};


//object 1 points
Point(9) = {0, 0, 0.5};
Point(10) = {1, 0, 0.5};
Point(11) = {0, 1, 0.5};
Point(12) = {1, 1, 0.5};

//object 1 lines
Line(13) = {9, 10};
Line(14) = {10, 12};
Line(15) = {12, 11};
Line(16) = {11, 9};

//object 1 surface
Curve Loop(8) = {13, 14, -15, -16};
Plane Surface(7) = {8};


//object 2 points
Point(13) = {0.5, 0, 0};
Point(14) = {0.5, 0, 1};
Point(15) = {0.5, 1, 1};
Point(16) = {0.5, 1, 0};

//object 2 lines
Line(17) = {13, 14};
Line(18) = {14, 15};
Line(19) = {15, 16};
Line(20) = {16, 13};

//object 2 surface
Curve Loop(9) = {17, 18, 19, -20};
Plane Surface(8) = {9};



//Center point definition
//Point(17) = {0.5, 0.5, 0.5, 1.0,h};
//Physical Point("Center_Point") = {17};

BooleanFragments{ Surface{7,8}; Delete; }{ }

//object 3 points
Point(21) = {0, 0.5, 0};
Point(22) = {0, 0.5, 1};
Point(23) = {1, 0.5, 1};
Point(24) = {1, 0.5, 0};

//object 3 lines
Line(26) = {21, 22};
Line(27) = {22, 23};
Line(28) = {23, 24};
Line(29) = {24, 21};

//object 3 surface
Curve Loop(12) = {26, 27, 28, -29};
Plane Surface(11) = {12};

bc_bac[]= {1,2,3,4,5,6};

BooleanFragments{ Surface{7,8,9,10,11,bc_bac[]}; Delete; }{ }


bc_top[] = {21,22,23,24};
bc_bottom[] = {15,13,14,16};
bc_front[] = {33,34,35,36};
bc_right[] = {17,18,19,20};
bc_left[] = {25,26,27,28};
bc_back[] = {29,30,31,32};
f1[] = {1,2,3,4};
f2[] = {5,6,7,8};
f3[] = {9,10,11,12};

bc[] = {bc_top[],bc_bottom[],bc_front[],bc_right[],bc_left[],bc_back[]};

//Domain volume
Surface Loop(1) = {bc[]};
Volume(1) = {1};

Surface {f1[],f2[],f3[],bc[]} In Volume {1};


Transfinite Line "*" = n+1;
Transfinite Surface "*";
Transfinite Volume "*";
Recombine Surface "*";
Recombine Volume "*";



Physical Volume("Domain_Vol") = {1};

Physical Surface("BC_Top") = {bc_top[]};
Physical Surface("BC_Bottom") = {bc_bottom[]};
Physical Surface("BC_Front") = {bc_front[]};
Physical Surface("BC_Right") = {bc_right[]};
Physical Surface("BC_Left") = {bc_left[]};
Physical Surface("BC_Back") = {bc_back[]};

Physical Surface("Frac_1") = {f1[]};
Physical Surface("Frac_2") = {f2[]};
Physical Surface("Frac_3") = {f3[]};


