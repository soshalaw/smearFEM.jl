//-------------------------------------------------------------
// Parameters (INPUT: modify these values)
//-------------------------------------------------------------
radius = 25.0;          // cylinder radius in mm
height = 40.0;          // cylinder height in mm
elem_size_2d = 10;      // characteristic element size

//-------------------------------------------------------------
// Bottom circle (z = 0)
//-------------------------------------------------------------
Point(1) = {0, 0, 0, elem_size_2d};
Point(2) = { radius, 0, 0, elem_size_2d};
Point(3) = { 0, radius, 0, elem_size_2d};
Point(4) = {-radius, 0, 0, elem_size_2d};
Point(5) = { 0,-radius, 0, elem_size_2d};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

//-------------------------------------------------------------
// Top circle (z = height)
//-------------------------------------------------------------
Point(6) = {0, 0, height, elem_size_2d};
Point(7) = { radius, 0, height, elem_size_2d};
Point(8) = { 0, radius, height, elem_size_2d};
Point(9) = {-radius, 0, height, elem_size_2d};
Point(10) = { 0,-radius, height, elem_size_2d};

Circle(5) = {7, 6, 8};
Circle(6) = {8, 6, 9};
Circle(7) = {9, 6, 10};
Circle(8) = {10, 6, 7};

Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};

//-------------------------------------------------------------
// Vertical lines connecting bottom to top rim
//-------------------------------------------------------------
Line(9)  = {2, 7};
Line(10) = {3, 8};
Line(11) = {4, 9};
Line(12) = {5, 10};

//-------------------------------------------------------------
// Lateral surfaces (four quarter-cylinder patches)
// Each patch is bounded by: bottom arc, vertical line, top arc, vertical line
//-------------------------------------------------------------
Curve Loop(3) = {1, 10, -5, -9};
Surface(3) = {3};

Curve Loop(4) = {2, 11, -6, -10};
Surface(4) = {4};

Curve Loop(5) = {3, 12, -7, -11};
Surface(5) = {5};

Curve Loop(6) = {4, 9, -8, -12};
Surface(6) = {6};

//-------------------------------------------------------------
// Volume
//-------------------------------------------------------------
Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

//-------------------------------------------------------------
// Physical groups
//-------------------------------------------------------------
Physical Surface("Bottom") = {1};
Physical Surface("Top") = {2};
Physical Surface("Lateral") = {3, 4, 5, 6};
Physical Volume("Cylinder") = {1};
