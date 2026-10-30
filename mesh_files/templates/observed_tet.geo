// Paper-2 initial mesh: volume inside the closed surface from smearFEM's `mesh_from_field`
// (surface.msh: entity 1 = cap on the table, 2 = free surface, shared nodes).
elem_size = 3.0;
Merge "surface.msh";
CreateTopology;
CreateGeometry;
Surface Loop(1) = {1, 2};
Volume(1) = {1};
Physical Surface("Bottom") = {1};
Physical Surface("Free") = {2};
Physical Volume("Volume") = {1};
Mesh.MeshSizeMax = elem_size;
Mesh.MeshSizeMin = elem_size / 4;
// mid-edge nodes on straight edges: curved onto the reparametrized surface they fold P2 tets
Mesh.SecondOrderLinear = 1;
