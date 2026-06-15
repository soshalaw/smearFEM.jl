//-------------------------------------------------------------
// Parameters (INPUT: modify these values)
//-------------------------------------------------------------
lx = 1.0;
ly = 1.0;
lz = 1.0;
elem_size = 0.25;

hx = lx / 2;
hy = ly / 2;
nx = Ceil(lx / elem_size);
ny = Ceil(ly / elem_size);
nz = Ceil(lz / elem_size);

//-------------------------------------------------------------
// Bottom face (z = 0) — centered at origin in XY
//-------------------------------------------------------------
Point(1) = {-hx, -hy, 0, elem_size};  // front-left
Point(2) = { hx, -hy, 0, elem_size};  // front-right
Point(3) = { hx,  hy, 0, elem_size};  // back-right
Point(4) = {-hx,  hy, 0, elem_size};  // back-left

Line(1) = {1, 2};  // Front (y=-hy)
Line(2) = {2, 3};  // Right (x=hx)
Line(3) = {3, 4};  // Back  (y=hy)
Line(4) = {4, 1};  // Left  (x=-hx)

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Transfinite Line{1, 3} = nx + 1;
Transfinite Line{2, 4} = ny + 1;
Transfinite Surface{1} = {1, 2, 3, 4};
Recombine Surface{1};

//-------------------------------------------------------------
// Extrude in z to create hex volume
//-------------------------------------------------------------
out[] = Extrude {0, 0, lz} {
  Surface{1};
  Layers{nz};
  Recombine;
};

//-------------------------------------------------------------
// Physical groups
// out[0]=Top  out[1]=Volume
// out[2]=Front(y=0)  out[3]=Right(x=lx)
// out[4]=Back(y=ly)  out[5]=Left(x=0)
//-------------------------------------------------------------
Physical Surface("Bottom") = {1};
Physical Surface("Top")    = {out[0]};
Physical Surface("Front")  = {out[2]};
Physical Surface("Right")  = {out[3]};
Physical Surface("Back")   = {out[4]};
Physical Surface("Left")   = {out[5]};
Physical Volume("Box")     = {out[1]};
