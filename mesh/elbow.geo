L = 1.5;
R = 1;
lc = 0.05;
ri = 0.5;
re = 1.5;
 
//+
Point(1) = {0, 0, 0, lc};
//+
Point(2) = {L, 0, 0, lc};
//+
Point(3) = {2 * L, L, 0, lc};
//+
Point(4) = {2 * L, 2 * L, 0, lc};
//+
Point(5) = {2 * L-R, 2 * L, 0, lc};
//+
Point(6) = {2 * L-R, L, 0, lc};
//+
Point(7) = {2 * L-R, L, 0, lc};
//+
Point(8) = {L, R, 0, lc};
//+
Point(9) = {0, R, 0, lc};
//+
Point(10) = {re, re, 0, lc};
//+
Line(1) = {1, 2};
//+
Line(2) = {3, 4};
//+
Line(3) = {4, 5};
//+
Line(4) = {5, 6};
//+
Line(5) = {8, 9};
//+
Line(6) = {9, 1};
//+
Circle(7) = {8, 10, 6};
//+
Circle(8) = {2, 10, 3};
//+
Physical Curve("inlet", 10) = {6};
//+
Physical Curve("walls", 20) = {4, 7, 5, 1, 8, 2};
//+
Physical Curve("outlet", 30) = {3};
//+
Curve Loop(1) = {5, 6, 1, 8, 2, 3, 4, -7};
//+
Plane Surface(1) = {1};
//+
Physical Surface("domain", 1) = {1};
//+
Transfinite Curve {6} = 50 Using Progression 1;
