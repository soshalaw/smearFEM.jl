//-------------------------------------------------------------
// Parameters (INPUT: modify these values)
//-------------------------------------------------------------
lx = 1.0;
ly = 1.0;
elem_size = 0.25;

hx = lx / 2;
hy = ly / 2;

//-------------------------------------------------------------
// 4 corners (XY plane, z=0) — centered at origin
//-------------------------------------------------------------
Point(1) = {-hx,  -hy, 0, elem_size};
Point(2) = { hx,  -hy, 0, elem_size};
Point(3) = { hx,   hy, 0, elem_size};
Point(4) = {-hx,   hy, 0, elem_size};

Line(1) = {1, 2};  // Bottom
Line(2) = {2, 3};  // Right
Line(3) = {3, 4};  // Top
Line(4) = {4, 1};  // Left

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

//-------------------------------------------------------------
// Physical groups
//-------------------------------------------------------------
Physical Curve("Bottom")   = {1};
Physical Curve("Right")    = {2};
Physical Curve("Top")      = {3};
Physical Curve("Left")     = {4};
Physical Surface("Square") = {1};
