//-------------------------------------------------------------
// Parameters (INPUT: modify these values)
//-------------------------------------------------------------
radius = 25.0;          // cylinder radius in mm
height = 40.0;          // cylinder height in mm
nz = 4;                 // number of elements along height

// 2D characteristic length: sized so n_2d ≈ nz² (gives ~cubic element distribution)
// Blossom on a Frontal-Delaunay circle mesh produces ~2.65× more quads than the
// ideal area formula (πR²/h²) due to curved-boundary refinement. Corrected:
//   n_2d ≈ 2.65·πR²/h² = nz²  →  h = R·√(2.65π)/nz
elem_size = radius * Sqrt(2.65 * Pi) / nz;

// Force pure-quad meshing of the base surface so extrusion gives pure hexes.
// SubdivisionAlgorithm is intentionally omitted: on a circular boundary Blossom
// leaves few unpaired triangles, but SubdivisionAlgorithm=1 would split each
// remaining triangle into 3 quads, inflating the element count ~6×.
Mesh.RecombineAll = 1;
Mesh.RecombinationAlgorithm = 2;

//-------------------------------------------------------------
// Points (center and rim)
//-------------------------------------------------------------
Point(0) = {0, 0, 0, elem_size};        // center
Point(1) = { radius, 0, 0, elem_size}; // right
Point(2) = { 0, radius, 0, elem_size}; // top
Point(3) = {-radius, 0, 0, elem_size}; // left
Point(4) = { 0, -radius, 0, elem_size}; // bottom

// Circle arcs
Circle(1) = {1, 0, 2};
Circle(2) = {2, 0, 3};
Circle(3) = {3, 0, 4};
Circle(4) = {4, 0, 1};

// Full circle surface
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Extrude 2D mesh to create 3D cylinder
out[] = Extrude {0, 0, height} {
  Surface{1};
  Layers{nz};
  Recombine;
};

// Physical groups
Physical Surface("Bottom") = {1};
Physical Surface("Top") = {out[0]};
Physical Surface("Lateral") = {out[2], out[3], out[4], out[5]};
Physical Volume("Volume") = {out[1]};