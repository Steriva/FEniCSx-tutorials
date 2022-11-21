lc = 1.;

scale = 1; // the geometry is made in cm

L = 38 / scale;

fuel_or = 1.791 / scale;
clad_ir = 1.800 / scale;
clad_or = 1.880 / scale;

R = 4.05 / 2 / scale; // pitch

// geometry: points
Point(1) = {-L/2, 0, 0, lc};
Point(2) = {L/2, 0, 0, lc};

Point(3) = {L/2, fuel_or, 0, lc};
Point(4) = {L/2, clad_ir, 0, lc};
Point(5) = {L/2, clad_or, 0, lc};

Point(6) = {L/2, R, 0, lc};
Point(7) = {-L/2, R, 0, lc};

Point(8)  = {-L/2, clad_or, 0, lc};
Point(9)  = {-L/2, clad_ir, 0, lc};
Point(10) = {-L/2, fuel_or, 0, lc};

// Lines
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 9};
//+
Line(9) = {9, 10};
//+
Line(10) = {10, 1};
//+
Line(11) = {10, 3};
//+
Line(12) = {4, 9};
//+
Line(13) = {8, 5};
//+
Curve Loop(1) = {1, 2, -11, 10};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {11, 3, 12, 9};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {8, -12, 4, -13};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {5, 6, 7, 13};
//+
Plane Surface(4) = {4};
//+
Physical Curve("out", 10) = {6};
//+
Physical Curve("sym", 20) = {1};
//+
Physical Curve("top", 30) = {5, 4, 3, 2};
//+
Physical Curve("bottom", 40) = {10, 9, 8};
//+
Physical Curve("inlet", 50) = {7};
//+
Physical Surface("fuel", 1) = {1};
//+
Physical Surface("helium", 2) = {2};
//+
Physical Surface("cladding", 3) = {3};
//+
Physical Surface("water", 4) = {4};
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
//+
Transfinite Surface {3};
//+
Transfinite Surface {4};
//+
Transfinite Curve {1, 12, 11, 6, 13} = 70 * 5 Using Progression 1;
//+
Transfinite Curve {10, 2} = 30 Using Progression 1;
//+
Transfinite Curve {9, 3} = 2 Using Progression 1;
//+
Transfinite Curve {4, 8} = 5 Using Progression 1;
//+
Transfinite Curve {7} = 20 Using Progression 0.925;
//+
Transfinite Curve {5} = 20 Using Progression 1/0.925;
//+
Recombine Surface {1, 3, 2, 4};
