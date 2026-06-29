//-------------------------------------------------------------
// Parameters (INPUT: modify these values)
//-------------------------------------------------------------
lx = 1.0;
ly = 1.0;
lz = 1.0;
r  = 0.1;         // vertical edge fillet radius
nz = 4;           // number of elements along height (primary input)

elem_size = lz / nz;

hx = lx / 2;
hy = ly / 2;

//-------------------------------------------------------------
// Bottom face points (z = 0) — rounded rectangle, centered at origin in XY
// 8 arc-endpoints + 4 arc centers, ordered CCW from above
//-------------------------------------------------------------
Point(1)  = {-hx+r,   -hy,    0, elem_size};  // front-left
Point(2)  = { hx-r,   -hy,    0, elem_size};  // front-right
Point(3)  = { hx,     -hy+r,  0, elem_size};  // right-front
Point(4)  = { hx,      hy-r,  0, elem_size};  // right-back
Point(5)  = { hx-r,    hy,    0, elem_size};  // back-right
Point(6)  = {-hx+r,    hy,    0, elem_size};  // back-left
Point(7)  = {-hx,      hy-r,  0, elem_size};  // left-back
Point(8)  = {-hx,     -hy+r,  0, elem_size};  // left-front
// arc centers (z = 0)
Point(9)  = {-hx+r,   -hy+r,  0, elem_size};
Point(10) = { hx-r,   -hy+r,  0, elem_size};
Point(11) = { hx-r,    hy-r,  0, elem_size};
Point(12) = {-hx+r,    hy-r,  0, elem_size};

//-------------------------------------------------------------
// Top face points (z = lz) — same layout
//-------------------------------------------------------------
Point(13) = {-hx+r,   -hy,    lz, elem_size};
Point(14) = { hx-r,   -hy,    lz, elem_size};
Point(15) = { hx,     -hy+r,  lz, elem_size};
Point(16) = { hx,      hy-r,  lz, elem_size};
Point(17) = { hx-r,    hy,    lz, elem_size};
Point(18) = {-hx+r,    hy,    lz, elem_size};
Point(19) = {-hx,      hy-r,  lz, elem_size};
Point(20) = {-hx,     -hy+r,  lz, elem_size};
// arc centers (z = lz)
Point(21) = {-hx+r,   -hy+r,  lz, elem_size};
Point(22) = { hx-r,   -hy+r,  lz, elem_size};
Point(23) = { hx-r,    hy-r,  lz, elem_size};
Point(24) = {-hx+r,    hy-r,  lz, elem_size};

//-------------------------------------------------------------
// Bottom edges
//-------------------------------------------------------------
Line(1)   = {1, 2};         // front  (y=0)
Circle(2) = {2, 10, 3};     // front-right arc
Line(3)   = {3, 4};         // right  (x=lx)
Circle(4) = {4, 11, 5};     // back-right arc
Line(5)   = {5, 6};         // back   (y=ly)
Circle(6) = {6, 12, 7};     // back-left arc
Line(7)   = {7, 8};         // left   (x=0)
Circle(8) = {8, 9,  1};     // front-left arc

//-------------------------------------------------------------
// Top edges
//-------------------------------------------------------------
Line(9)    = {13, 14};
Circle(10) = {14, 22, 15};
Line(11)   = {15, 16};
Circle(12) = {16, 23, 17};
Line(13)   = {17, 18};
Circle(14) = {18, 24, 19};
Line(15)   = {19, 20};
Circle(16) = {20, 21, 13};

//-------------------------------------------------------------
// Vertical edges (one per arc endpoint)
//-------------------------------------------------------------
Line(17) = {1,  13};
Line(18) = {2,  14};
Line(19) = {3,  15};
Line(20) = {4,  16};
Line(21) = {5,  17};
Line(22) = {6,  18};
Line(23) = {7,  19};
Line(24) = {8,  20};

//-------------------------------------------------------------
// Surfaces
//-------------------------------------------------------------
// Bottom cap
Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
Plane Surface(1) = {1};

// Top cap
Curve Loop(2) = {9, 10, 11, 12, 13, 14, 15, 16};
Plane Surface(2) = {2};

// Flat side faces
Curve Loop(3) = {1,  18, -9,  -17};   Plane Surface(3) = {3};  // Front  y=0
Curve Loop(4) = {3,  20, -11, -19};   Plane Surface(4) = {4};  // Right  x=lx
Curve Loop(5) = {5,  22, -13, -21};   Plane Surface(5) = {5};  // Back   y=ly
Curve Loop(6) = {7,  24, -15, -23};   Plane Surface(6) = {6};  // Left   x=0

// Quarter-cylinder corner surfaces (curved — uses surface filling)
Curve Loop(7)  = {2,  19, -10, -18};  Surface(7)  = {7};   // Front-right
Curve Loop(8)  = {4,  21, -12, -20};  Surface(8)  = {8};   // Back-right
Curve Loop(9)  = {6,  23, -14, -22};  Surface(9)  = {9};   // Back-left
Curve Loop(10) = {8,  17, -16, -24};  Surface(10) = {10};  // Front-left

Surface Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
Volume(1) = {1};

//-------------------------------------------------------------
// Physical groups
//-------------------------------------------------------------
Physical Surface("Bottom") = {1};
Physical Surface("Top")    = {2};
Physical Surface("Front")  = {3, 7};   // flat front + front-right arc
Physical Surface("Right")  = {4, 8};   // flat right + back-right arc
Physical Surface("Back")   = {5, 9};   // flat back + back-left arc
Physical Surface("Left")   = {6, 10};  // flat left + front-left arc
Physical Volume("Box")     = {1};
