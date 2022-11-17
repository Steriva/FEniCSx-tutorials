SetFactory("OpenCASCADE");

lc = 1;
R = 0.5;
L = 10;


Point(1) = {0, 0, 0, lc};
//+
Point(2) = {L, -R, 0, lc};
//+
Point(3) = {L, R, 0, lc};
//+
Point(4) = {0, R, 0, lc};
//+
Point(5) = {0, -R, 0, lc};
//+
Line(1) = {5, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};

Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Physical Curve("inlet", 10) = {4};
//+
Physical Curve("top", 20) = {3};
//+
Physical Curve("bottom", 30) = {1};
//+
Physical Curve("outlet", 40) = {2};
//+
Physical Surface("domain", 1) = {1};

n = 10;
//+
Transfinite Curve {1, 3} = 12*n Using Progression 1;
//+
Transfinite Curve {2} = 4*n Using Progression 1;
//+
Transfinite Curve {4} = 4*n Using Progression 1;
//+
Transfinite Surface {1};
//+
Recombine Surface {1};
