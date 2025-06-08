// Gmsh project created on Sun Jun  8 16:43:53 2025
SetFactory("OpenCASCADE");
//+
Point(1) = {1.0, 0.0, 0, 1.0};
//+
Point(2) = {1.0, 1.0, 0, 1.0};
//+
Point(3) = {1.0, 1.0, 1.0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Point(5) = {0, 0, 0, 1.0};
//+
Point(6) = {0, 0, 1, 1.0};
//+
Point(7) = {1, 0, 1, 1.0};
//+
Point(8) = {0, 1, 1, 1.0};
//+
Line(1) = {5, 1};
//+
Line(2) = {1, 2};
//+
Line(3) = {2, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {6, 8};
//+
Line(6) = {8, 3};
//+
Line(7) = {3, 7};
//+
Line(8) = {7, 6};
//+
Line(9) = {6, 5};
//+
Line(10) = {7, 1};
//+
Line(11) = {3, 2};
//+
Line(12) = {8, 4};
//+
Curve Loop(1) = {7, 10, 2, -11};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {6, 11, 3, -12};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {5, 12, 4, -9};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {9, 1, -10, 8};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {4, 1, 2, 3};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {6, 7, 8, 5};
//+
Plane Surface(6) = {6};
//+
Surface Loop(7) = {4, 3, 6, 2, 1, 5};
//+
Volume(2) = {7};
