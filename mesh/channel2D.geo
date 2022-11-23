lc = 1;
scale = 1; // the geometry is made in cm

L_active = 35.56 / scale;
L_graphi = 10.2  / scale;
L_clad_o = 1.27  / scale;

fuel_or = 1.791 / scale;
clad_ir = 1.800 / scale;
clad_or = 1.880 / scale;

L = L_active + 2 * L_graphi + 2 * L_clad_o;

R = 4.36 / 2 / scale; // pitch

// geometry: points

Point(5) = {L/2, clad_or, 0, lc};

Point(6) = {L/2, R, 0, lc};
Point(7) = {-L/2, R, 0, lc};

Point(8)  = {-L/2, clad_or, 0, lc};

// Lines
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 8};
//+
Line(13) = {8, 5};
//+
Physical Curve("inlet", 10) = {7};
//+
Physical Curve("bottom", 20) = {13};
//+
Physical Curve("top", 30) = {6};
//+
Physical Curve("out", 40) = {5};
//+
Curve Loop(1) = {6, 7, 13, 5};
//+
Plane Surface(1) = {1};
//+
Physical Surface("domain", 1) = {1};

n = 10;  
//+
Transfinite Curve {13, 6} = 30*n Using Progression 1;
//+
Transfinite Curve {7} = 2*n Using Progression 0.95;
//+
Transfinite Curve {5} = 2*n Using Progression 1/0.95;
//+
Transfinite Surface {1};
//+
Recombine Surface {1};