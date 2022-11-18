L = 1.5;
R = 1;
lc = 0.036;
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
Point(11) = {2, 1, 0, lc};
//+
Point(12) = {3, 0, 0, lc};
//+
Line(7) = {2, 12};
//+
Line(8) = {12, 3};
//+
Line(9) = {8, 11};
//+
Line(10) = {11, 6};
//+
Curve Loop(1) = {5, 6, 1, 7, 8, 2, 3, 4, -10, -9};
//+
Plane Surface(1) = {1};
//+
Physical Curve("inlet", 10) = {6};
//+
Physical Curve("walls", 20) = {7, 1, 5, 9, 10, 4, 2, 8};
//+
Physical Curve("outlet", 30) = {3};
//+
Physical Surface("domain", 1) = {1};
