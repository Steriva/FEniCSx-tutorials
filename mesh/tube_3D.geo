lc = 3.5e-2; // the characteristic length of the mesh (to be expanded by m4)

bound = 0.5;

// geometry
Point(1) = {0, 0, 0, lc};
Point(2) = {bound, 0, 0, lc};
Point(3) = {bound, bound, 0, lc};
Point(4) = {0, bound, 0, lc};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Extrude {0, 0, 2} {
  Surface{1}; 
}
//+
Physical Surface("inlet", 10) = {1};
//+
Physical Surface("walls", 20) = {25, 21, 13, 17};
//+
Physical Surface("outlet", 30) = {26};
//+
Physical Volume("domain", 1) = {1};
