//-------------------------------------------------------------
// Parameters (INPUT: modify these values)
//-------------------------------------------------------------
length = 1.0;
elem_size = 0.25;

n = Ceil(length / elem_size) + 1;

//-------------------------------------------------------------
// Two endpoints along x-axis — centered at origin
//-------------------------------------------------------------
Point(1) = {-length/2, 0, 0, elem_size};
Point(2) = { length/2, 0, 0, elem_size};

Line(1) = {1, 2};
Transfinite Line {1} = n;

//-------------------------------------------------------------
// Physical groups
//-------------------------------------------------------------
Physical Point("Left")  = {1};
Physical Point("Right") = {2};
Physical Curve("Line")  = {1};
