//-------------------------------------------------------------
// Parameters (INPUT: modify these values)
//-------------------------------------------------------------
lx = 1.0;
ly = 1.0;
lz = 1.0;
r  = 0.1;         // vertical edge fillet radius
nz = 4;           // number of elements along height (primary input)

hx = lx / 2;
hy = ly / 2;
// Full-quad recombination (algorithm 2) inflates n_2D by ~2x on a rounded-rectangle
// base (not the ~3.5x assumed previously, which left the mesh half as fine as asked).
// Multiply char length by sqrt(2) so the post-subdivision count stays near nz^2,
// giving ~nz^3 elements after extrusion, matching the cylinder template's convention.
// The recombiner quantises n_2D onto plateaus (48 / 124 / 272), so this factor sits
// mid-plateau rather than at a boundary.
// NOTE: do not write "<param> = value" inside these comments. _generate_mesh_geo
// substitutes params with the regex `<param>\s*=\s*[^;]+;`, which spans newlines and
// would swallow everything up to the next semicolon, including this definition.
elem_size_xy = Sqrt(2.0 * lx * ly) / nz;

//-------------------------------------------------------------
// Bottom face (z = 0) — rounded rectangle, centered at origin in XY
// 8 arc-endpoints + 4 arc centers, ordered CCW from above
//-------------------------------------------------------------
Point(1)  = {-hx+r,   -hy,    0, elem_size_xy};  // front-left
Point(2)  = { hx-r,   -hy,    0, elem_size_xy};  // front-right
Point(3)  = { hx,     -hy+r,  0, elem_size_xy};  // right-front
Point(4)  = { hx,      hy-r,  0, elem_size_xy};  // right-back
Point(5)  = { hx-r,    hy,    0, elem_size_xy};  // back-right
Point(6)  = {-hx+r,    hy,    0, elem_size_xy};  // back-left
Point(7)  = {-hx,      hy-r,  0, elem_size_xy};  // left-back
Point(8)  = {-hx,     -hy+r,  0, elem_size_xy};  // left-front
// arc centers
Point(9)  = {-hx+r,   -hy+r,  0, elem_size_xy};
Point(10) = { hx-r,   -hy+r,  0, elem_size_xy};
Point(11) = { hx-r,    hy-r,  0, elem_size_xy};
Point(12) = {-hx+r,    hy-r,  0, elem_size_xy};

Line(1)   = {1, 2};
Circle(2) = {2, 10, 3};
Line(3)   = {3, 4};
Circle(4) = {4, 11, 5};
Line(5)   = {5, 6};
Circle(6) = {6, 12, 7};
Line(7)   = {7, 8};
Circle(8) = {8, 9,  1};

Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8};
Plane Surface(1) = {1};
Recombine Surface {1};

Mesh.RecombineAll = 1;
Mesh.RecombinationAlgorithm = 2;

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
// out[2..9] = lateral surfaces in curve-loop order:
//   [2] Front(y=0)  [3] Front-right arc
//   [4] Right(x=lx) [5] Back-right arc
//   [6] Back(y=ly)  [7] Back-left arc
//   [8] Left(x=0)   [9] Front-left arc
//-------------------------------------------------------------
Physical Surface("Bottom") = {1};
Physical Surface("Top")    = {out[0]};
Physical Surface("Front")  = {out[2], out[3]};  // flat front + front-right arc
Physical Surface("Right")  = {out[4], out[5]};  // flat right + back-right arc
Physical Surface("Back")   = {out[6], out[7]};  // flat back + back-left arc
Physical Surface("Left")   = {out[8], out[9]};  // flat left + front-left arc
Physical Volume("Box")     = {out[1]};
