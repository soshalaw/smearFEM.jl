//-------------------------------------------------------------
// Parameters (INPUT: modify these values)
//-------------------------------------------------------------
radius = 25.0;          // cylinder radius in mm
height = 40.0;          // cylinder height in mm
elem_size_2d = 10;      // characteristic element size (2D plane)

// Calculated parameter
nz = Ceil(height / elem_size_2d);  // number of layers in z-direction (auto-calculated)

//-------------------------------------------------------------
// Points (center and rim)
//-------------------------------------------------------------
Point(0) = {0, 0, 0, elem_size_2d};        // center
Point(1) = { radius, 0, 0, elem_size_2d}; // right
Point(2) = { 0, radius, 0, elem_size_2d}; // top
Point(3) = {-radius, 0, 0, elem_size_2d}; // left
Point(4) = { 0, -radius, 0, elem_size_2d}; // bottom

// Circle arcs
Circle(1) = {1, 0, 2};
Circle(2) = {2, 0, 3};
Circle(3) = {3, 0, 4};
Circle(4) = {4, 0, 1};

// Full circle surface
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Recombine to create quad mesh (unstructured)
Recombine Surface {1};

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