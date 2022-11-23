lc = 0.1;

scale = 1; // the geometry is made in cm

L_active = 35.56 / scale;
L_graphi = 10.2  / scale;
L_clad_o = 1.27  / scale;

fuel_or = 1.791 / scale;
clad_ir = 1.800 / scale;
clad_or = 1.880 / scale;

R = 4.36 / 2 / scale; // pitch

// geometry: points
Point(1) = {-L_active/2 - L_graphi - L_clad_o, 0, 0, lc};
Point(2) = {-L_active/2 - L_graphi           , 0, 0, lc};
Point(3) = {-L_active/2                      , 0, 0, lc};

Point(4) = {L_active/2 + L_graphi + L_clad_o, 0, 0, lc};
Point(5) = {L_active/2 + L_graphi           , 0, 0, lc};
Point(6) = {L_active/2                      , 0, 0, lc};

Point(7) = {-L_active/2 - L_graphi - L_clad_o, fuel_or, 0, lc};
Point(8) = {-L_active/2 - L_graphi           , fuel_or, 0, lc};
Point(9) = {-L_active/2                      , fuel_or, 0, lc};

Point(10) = {L_active/2 + L_graphi + L_clad_o, fuel_or, 0, lc};
Point(11) = {L_active/2 + L_graphi           , fuel_or, 0, lc};
Point(12) = {L_active/2                      , fuel_or, 0, lc};

Point(13) = {-L_active/2 - L_graphi - L_clad_o, clad_ir, 0, lc};
Point(14) = { L_active/2 + L_graphi + L_clad_o, clad_ir, 0, lc};

Point(15) = {-L_active/2 - L_graphi - L_clad_o, clad_or, 0, lc};
Point(16) = { L_active/2 + L_graphi + L_clad_o, clad_or, 0, lc};

Point(17) = {-L_active/2 - L_graphi - L_clad_o, R, 0, lc};
Point(18) = { L_active/2 + L_graphi + L_clad_o, R, 0, lc};

// Lines

Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 6};
//+
Line(4) = {6, 5};
//+
Line(5) = {5, 4};
//+
Line(6) = {7, 8};
//+
Line(7) = {8, 9};
//+
Line(8) = {9, 12};
//+
Line(9) = {12, 11};
//+
Line(10) = {11, 10};
//+
Line(11) = {13, 14};
//+
Line(12) = {15, 16};
//+
Line(13) = {17, 18};
//+
Line(14) = {1, 7};
//+
Line(15) = {2, 8};
//+
Line(16) = {3, 9};
//+
Line(17) = {6, 12};
//+
Line(18) = {5, 11};
//+
Line(19) = {4, 10};
//+
Line(20) = {7, 13};
//+
Line(21) = {10, 14};
//+
Line(22) = {13, 15};
//+
Line(23) = {14, 16};
//+
Line(24) = {15, 17};
//+
Line(25) = {16, 18};

// Surfaces
Curve Loop(1) = {1, 15, -6, -14};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {5, 19, -10, -18};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {11, 23, -12, -22};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {7, -16, -2, 15};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {4, 18, -9, -17};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {6, 7, 8, 9, 10, 21, -11, -20};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {3, 17, -8, -16};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {13, -25, -12, 24};
//+
Plane Surface(8) = {8};

// Physical Groups

Physical Surface("fuel", 1) = {7};
//+
Physical Surface("helium", 2) = {6};
//+
Physical Surface("cladding", 3) = {3, 1, 2};
//+
Physical Surface("water", 4) = {8};
//+
Physical Surface("graphite", 5) = {4, 5};
//+
Physical Curve("middle", 10) = {13};
//+
Physical Curve("outlet", 20) = {19, 21, 23, 25};
//+
Physical Curve("centre", 30) = {5, 4, 3, 1, 2};
//+
Physical Curve("inlet", 40) = {14, 20, 22, 24};

// Mesh
Recombine Surface {1, 2, 3, 4, 5, 6, 7, 8};
//+
Transfinite Surface {1,2,3,4,5,7,8};
//+
Transfinite Curve {24} = 15 Using Progression 1/0.85;
//+
Transfinite Curve {25} = 15 Using Progression 0.85;