SetFactory("OpenCASCADE");

Point(1)  = { -1,  -1, 0, 1.0};
Point(2)  = {  1,  -1, 0, 1.0};
Point(3)  = {  1,   1, 0, 1.0};
Point(4)  = { -1,   1, 0, 1.0};

Mesh.MeshSizeFactor = 0.1;
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Circle(5) = {0., 0.0, 0, 0.5, 0, 2*Pi};

//+
Curve Loop(1) = {5};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {4, 1, 2, 3};
//+
Curve Loop(3) = {5};
//+
Plane Surface(2) = {2, 3};
//+
Physical Surface("omega_1", 1) = {1};
//+
Physical Surface("omega_2", 2) = {2};
//+
Physical Curve("top", 10) = {3};
//+
Physical Curve("side", 20) = {2, 4};
//+
Physical Curve("bottom", 30) = {1};
