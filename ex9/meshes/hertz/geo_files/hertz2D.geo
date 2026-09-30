// -----------------------------------------------------------------------------
// Parameters & Geometry for 2D Hertz Contact (Cylinder on Flat)
// -----------------------------------------------------------------------------
R  = 10.0; // Indenter radius
W  = 4.0;  // Half-width (using symmetry)
Hf = 4.0;  // Foundation depth
Hi = 4.0;  // Indenter height
lc = 0.2;  // Mesh characteristic length

// Foundation
Point(1) = {0, -Hf, 0, lc};
Point(2) = {W, -Hf, 0, lc};
Point(3) = {W, 0, 0, lc};
Point(4) = {0, 0, 0, lc}; // Contact point

Line(1) = {1, 2}; // FOUND_BOT
Line(2) = {2, 3}; // FOUND_RIGHT
Line(3) = {3, 4}; // FOUND_TOP
Line(4) = {4, 1}; // SYM_FOUND (x=0)

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1}; 

// Void (Third Medium Gap)
y_gap = R - Sqrt(R*R - W*W);
Point(5) = {W, y_gap, 0, lc};
Point(100) = {0, R, 0, lc}; // Sphere center
Circle(5) = {5, 100, 4};    // IND_BOT (Arc)

Line(6) = {3, 5}; // VOID_RIGHT
Curve Loop(2) = {6, 5, -3}; // Correctly oriented CCW
Plane Surface(2) = {2};

// Indenter
Point(6) = {W, Hi, 0, lc};
Point(7) = {0, Hi, 0, lc};
Line(7) = {5, 6}; // IND_RIGHT
Line(8) = {6, 7}; // IND_TOP
Line(9) = {7, 4}; // SYM_IND (x=0)

Curve Loop(3) = {7, 8, 9, -5};
Plane Surface(3) = {3}; 

Coherence;

// -----------------------------------------------------------------------------
// Physical Groups Output
// -----------------------------------------------------------------------------
Physical Surface("FOUNDATION", 100) = {1};
Physical Surface("VOID", 101) = {2};
Physical Surface("INDENTER", 102) = {3};

Physical Curve("FOUND_BOT", 200) = {1};
Physical Curve("SYM_FOUND", 201) = {4};
Physical Curve("SYM_IND", 202) = {9};
Physical Curve("IND_TOP", 203) = {8};

// -----------------------------------------------------------------------------
// Quadrilateral Meshing Algorithms[cite: 1]
// -----------------------------------------------------------------------------
Transfinite Curve {1, 3, 5, 8} = 40;
Transfinite Curve {2, 4} = 20;
Transfinite Curve {6} = 5; // High resolution through the thin gap
Transfinite Curve {7, 9} = 20;

Transfinite Surface {1}; 
Transfinite Surface {2};
Transfinite Surface {3};
Recombine Surface {1, 2, 3}; 

Mesh.Algorithm = 8; 
Mesh.RecombinationAlgorithm = 1;
Mesh.Smoothing = 100;