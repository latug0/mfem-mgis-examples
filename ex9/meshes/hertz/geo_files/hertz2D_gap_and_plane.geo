// -----------------------------------------------------------------------------
// Parameters & Geometry for 2D Hertz Contact (Cylinder on Flat)
// -----------------------------------------------------------------------------
R  = 1; // Indenter radius
W  = 3*R;  // Half-width of plane (using symmetry)
// W >= R needed or >= (Pi/2)*R ?
W = Max(W,R);
Hf = 0.9*R;  // Foundation/plane depth
void_gap = 0.05*R; // Gap between indenter and foundation
arc_ratio = 0.9; // Percentage of arc length covered by third Medium, ie. contact surface on circle/indenter (between 0 and 1)
//plane_ratio = Sin(arc_ratio*Pi/2)*(R/W) ; // Percentage of plane length covered by third medium, ie.contact surface on plane/foundation

lc = 0.05;  // Mesh characteristic length





// Foundation
Point(1) = {0, -Hf, 0, lc};
Point(2) = {W, -Hf, 0, lc};
Point(3) = {W, 0, 0, lc};
Point(4) = {R*Sin(arc_ratio*Pi/2),0,0,lc}; //Same point as the third_medium limit on the circle arc
Point(5) = {0, 0, 0, lc}; // Contact point

Line(1) = {1, 2}; // FOUND_BOT
Line(2) = {2, 3}; // FOUND_RIGHT
Line(3) = {3, 4}; // FOUND_TOP
Line(4) = {4,  5};// FOUND_CONTACT
Line(5) = {5, 1}; // SYM_FOUND (x=0)

Curve Loop(1) = {1, 2, 3, 4, 5};
Plane Surface(1) = {1}; 

// Void (Third Medium Gap)

Point(6) = {0,void_gap,0,lc};
// Intermediate point to delimit contact surface on circle/indenter
Point(7) = {R*Sin(arc_ratio*Pi/2) , void_gap + R*(1-Cos(arc_ratio*Pi/2)), 0, lc};
Point(8) = {R, void_gap+R, 0, lc}; // Half_sphere
Point(9) = {0, R+void_gap, 0, lc}; // Sphere center

Circle(6) = {6, 9, 7};    // IND_CONTACT (Arc) // Zone covered by third medium
Circle(7) = {7, 9, 8};    // IND_BOT (Arc) // Zone not covered by third medium


Line(8) = {5,6}; // VOID_LEFT 
Line(9) = {7, 4}; // VOID_RIGHT
Curve Loop(2) = {8, 6, 9, 4}; // Third medium zone for the contact
Plane Surface(2) = {2};

// Indenter
Line(10) = {8, 9}; // IND_TOP
Line(11) = {9, 6}; // SYM_IND (x=0)

Curve Loop(3) = {6, 7, 10, 11}; // Circle/indenter zone
Plane Surface(3) = {3}; 

Coherence;

// -----------------------------------------------------------------------------
// Physical Groups Output
// -----------------------------------------------------------------------------
Physical Surface("FOUNDATION", 100) = {1};
Physical Surface("VOID", 101) = {2};
Physical Surface("INDENTER", 102) = {3};


Physical Curve("FOUND_BOT", 200) = {1};
Physical Curve("SYM_FOUND", 201) = {5};
Physical Curve("SYM_IND", 202) = {11};
Physical Curve("IND_TOP", 203) = {10};


//Line(1) = {1, 2}; // FOUND_BOT
//Line(2) = {2, 3}; // FOUND_RIGHT
//Line(3) = {3, 4}; // FOUND_TOP
//Line(4) = {4,  5};// FOUND_CONTACT
//Line(5) = {5, 1}; // SYM_FOUND (x=0)
//Circle(6) = {6, 9, 7};    // IND_CONTACT (Arc) // Zone covered by third medium
//Circle(7) = {7, 9, 8};    // IND_BOT (Arc) // Zone not covered by third medium
//Line(8) = {5,6}; // VOID_LEFT 
//Line(9) = {7, 4}; // VOID_RIGHT
//Line(10) = {8, 9}; // IND_TOP
//Line(11) = {6, 9}; // SYM_IND (x=0)
// -----------------------------------------------------------------------------
// Quadrilateral Meshing Algorithms[cite: 1]
// -----------------------------------------------------------------------------
//Transfinite Curve {1, 3, 5, 8} = 40;
//Transfinite Curve {2, 4} = 20;
//Transfinite Curve {6} = 5; // High resolution through the thin gap
//Transfinite Curve {7, 9} = 20;
Transfinite Curve {4,6} = 40; // Number of nodes in contact area
//Transfinite Curve {6} = Floor(40*Pi/2); // Number of nodes in contact area (scaled by quarter circle length)

Transfinite Curve {8,9} = 2; // 1 element width on for void
//Transfinite Surface {1}; 
Transfinite Surface {2}; // Void surface
Recombine Surface {2};
//Transfinite Surface {3};
//Recombine Surface {1, 2, 3}; 

//Mesh.Algorithm = 8; 
//Mesh.RecombinationAlgorithm = 1;
//Mesh.Smoothing = 100;