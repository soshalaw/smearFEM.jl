//-------------------------------------------------------------
// Parameters (INPUT: modify these values)
//-------------------------------------------------------------
radius = 1.0;
elem_size = 0.25;

//-------------------------------------------------------------
// Center + 4 rim points (XY plane, z=0)
// Quad mesh via Recombine on an extruded cross pattern
//-------------------------------------------------------------
Point(1) = { 0,       0,      0, elem_size};  // center
Point(2) = { radius,  0,      0, elem_size};
Point(3) = { 0,       radius, 0, elem_size};
Point(4) = {-radius,  0,      0, elem_size};
Point(5) = { 0,      -radius, 0, elem_size};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Recombine to produce quad elements
Recombine Surface {1};

//-------------------------------------------------------------
// Physical groups
//-------------------------------------------------------------
Physical Curve("Boundary") = {1, 2, 3, 4};
Physical Surface("Disk")   = {1};
