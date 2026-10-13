// -----------------------------------------------------------------------------
// Parameters & Geometry for 3D Hertz Contact (Sphere on Flat)
// -----------------------------------------------------------------------------
SetFactory("OpenCASCADE");
R  = 10.0;
W  = 3*4.0; 
Hf = 4.0;
Hi = 4.0;
lc = 0.5;

// Foundation
v_found = newv; Box(v_found) = {0, 0, -Hf, W, W, Hf};

// Bounding box for the upper region and the full Sphere
v_box = newv; Box(v_box) = {0, 0, 0, W, W, Hi};
v_sph = newv; Sphere(v_sph) = {0, 0, R, R};

// Indenter is the intersection of the Box and the Sphere
v_ind = newv; BooleanIntersection(v_ind) = { Volume{v_box}; Delete; }{ Volume{v_sph}; Delete; };

// Void is the remainder of the bounding box
v_box2 = newv; Box(v_box2) = {0, 0, 0, W, W, Hi};
v_void = newv; BooleanDifference(v_void) = { Volume{v_box2}; Delete; }{ Volume{v_ind}; };

// Fragment to ensure matching meshes at all boundaries
v_all() = BooleanFragments{ Volume{v_found, v_ind, v_void}; Delete; }{};

// -----------------------------------------------------------------------------
// 3D Physical Groups Output
// -----------------------------------------------------------------------------
Physical Volume("FOUNDATION") = {v_all(0)};
Physical Volume("INDENTER") = {v_all(1)};
Physical Volume("VOID") = {v_all(2)};

eps = 1e-3;
// Automatically group boundaries using bounding box limits
surf_found_bot() = Surface In BoundingBox {-eps, -eps, -Hf-eps, W+eps, W+eps, -Hf+eps};
Physical Surface("FOUND_BOT") = surf_found_bot[];

surf_sym_x() = Surface In BoundingBox {-eps, -eps, -Hf-eps, eps, W+eps, Hi+eps};
Physical Surface("SYM_X") = surf_sym_x[];

surf_sym_y() = Surface In BoundingBox {-eps, -eps, -Hf-eps, W+eps, eps, Hi+eps};
Physical Surface("SYM_Y") = surf_sym_y[];

surf_ind_top() = Surface In BoundingBox {-eps, -eps, Hi-eps, W+eps, W+eps, Hi+eps};
Physical Surface("IND_TOP") = surf_ind_top[];

Mesh.CharacteristicLengthMax = lc;
Mesh.CharacteristicLengthMin = lc/5;
Mesh.Algorithm3D = 10;
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;