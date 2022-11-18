SetFactory("OpenCASCADE");

lc = 1;
L = 2.2;
R_up = 0.16;
R_down = 0.15;
D = 0.15;
Radius = 0.05;

Point(1) = {0, 0, 0, lc};
//+
Point(2) = {L, -Radius-R_down, 0, lc};
//+
Point(3) = {L, Radius+R_up, 0, lc};
//+
Point(4) = {-D, Radius+R_up, 0, lc};
//+
Point(5) = {-D, -Radius-R_down, 0, lc};
//+
Line(1) = {5, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Circle(5) = {0, 0, 0, Radius, 0, 2*Pi};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Curve Loop(2) = {5};
//+
Plane Surface(1) = {1, 2};
//+
Physical Curve("inlet", 10) = {4};
//+
Physical Curve("wall", 20) = {5};
//+
Physical Curve("sym", 30) = {3, 1};
//+
Physical Curve("outlet", 40) = {2};
//+
Physical Surface("domain", 1) = {1};
//+
// Recombine Surface {1};
n = 20;
//+
Transfinite Curve {1} = 6*n Using Progression 1;
//+
Transfinite Curve {3} = 6*n Using Progression 1;
//+
Transfinite Curve {2} = n Using Progression 1;
//+
Transfinite Curve {4} = n Using Progression 1;
//+
Transfinite Curve {5} = 3*n Using Progression 1;
