// Gmsh project created on Fri Mar 22 17:07:30 2019
SetFactory("OpenCASCADE");

h=0.9;
Geometry.Tolerance = 1e-05;
Mesh.Algorithm = 8;

//Back side domain points
Point(1)={0,0,0,h};
Point(2)={1,0,0,h};
Point(3)={1,1,0,h};
Point(4)={0,1,0,h};

//Back side domain lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

//Back side domain surfaces
Curve Loop(1) = {3, 4, 1, 2};
Plane Surface(1) = {1};
Physical Surface("BC_Back") = {1};

//Right side domain points
Point(5)={1,0,1,h};
Point(6)={1,1,1,h};

//Right side domain lines
Line(5) = {2, 5};
Line(6) = {5, 6};
Line(7) = {6, 3};

//Right side domain surfaces
Curve Loop(2) = {2, -7, -6, -5};
Plane Surface(2) = {2};
Physical Surface("BC_Right") = {2};

//Front side domain points
Point(7)={0,0,1,h};
Point(8)={0,1,1,h};

//Front side domain lines
Line(8) = {6, 8};
Line(9) = {8, 7};
Line(10) = {7, 5};

//Front side domain surfaces
Curve Loop(3) = {6, 8, -9, -10};
Plane Surface(3) = {3};
Physical Surface("BC_Front") = {3};

//Left side domain points (NONE)

//Left side domain lines
Line(11) = {8, 4};
Line(12) = {1, 7};

//Left side domain surfaces
Curve Loop(4) = {-9, -11, -4, 12};
Plane Surface(4) = {4};
Physical Surface("BC_Left") = {4};

//Top side domain points (NONE)

//Top side domain lines (NONE)

//Top side domain surfaces
Curve Loop(6) = {3, -11, -8, 7};
Plane Surface(5) = {6};
Physical Surface("BC_Top") = {5};

//Bottom side domain points (NONE)

//Bottom side domain lines (NONE)

//Bottom side domain surfaces
Curve Loop(7) = {1, 5, -10, -12};
Plane Surface(6) = {7};
Physical Surface("BC_Bottom") = {6};

//Domain volume
Surface Loop(1) = {4, 3, 2, 1, 5, 6};
Volume(1) = {1};
Physical Volume("Domain_Vol") = {1};

//Fracture 1 points
Point(9) = {0, 0, 0.5, 1.0};
Point(10) = {1, 0, 0.5, 1.0};
Point(11) = {0, 1, 0.5, 1.0};
Point(12) = {1, 1, 0.5, 1.0};

//Fracture 1 lines
Line(13) = {9, 10};
Line(14) = {10, 12};
Line(15) = {12, 11};
Line(16) = {11, 9};

//Fracture 1 surface
Curve Loop(8) = {13, 14, -15, -16};
Plane Surface(7) = {8};
Physical Surface("Frac_1") = {7};

//Fracture 2 points
Point(13) = {0.5, 0, 0, 1.0};
Point(14) = {0.5, 0, 1, 1.0};
Point(15) = {0.5, 1, 1, 1.0};
Point(16) = {0.5, 1, 0, 1.0};

//Fracture 2 lines
Line(17) = {13, 14};
Line(18) = {14, 15};
Line(19) = {15, 16};
Line(20) = {16, 13};

//Fracture 2 surface
Curve Loop(9) = {17, 18, 19, -20};
Plane Surface(8) = {9};
Physical Surface("Frac_2") = {8};

//Center point definition
Point(17) = {0.5, 0.5, 0.5, 1.0,h};
Physical Point("Center_Point") = {17};
