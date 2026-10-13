// -----------------------------------------------------------------------------
// Parameters 
// -----------------------------------------------------------------------------
L = 1.0;
t = 0.072;
H = L / 2.0;

// Mesh Density
If (!Exists(ref))
  ref = 1;
EndIf
lc = (0.02 / ref) * L;

// -----------------------------------------------------------------------------
// Points Construction
// -----------------------------------------------------------------------------
// Solid C-Shape Domain Points
Point(1) = {0, 0, 0, lc};
Point(2) = {L, 0, 0, lc};
Point(3) = {L, t, 0, lc};
Point(4) = {t, t, 0, lc};
Point(5) = {t, H-t, 0, lc};
Point(6) = {L, H-t, 0, lc};
Point(7) = {L, H, 0, lc};
Point(8) = {L-0.1*t, H, 0, lc};
Point(9) = {0, H, 0, lc};

// Void Domain Expansion Points (spanning full height H at the right)
Point(10) = {L+1*t, 0, 0, lc};
Point(11) = {L+1*t, H, 0, lc};

// -----------------------------------------------------------------------------
// Curves Construction
// -----------------------------------------------------------------------------
Line(1) = {1, 2}; 
Line(2) = {2, 3}; // Right edge of bottom beam
Line(3) = {3, 4}; // Inner bottom edge
Line(4) = {4, 5}; // Inner vertical spine
Line(5) = {5, 6}; // Inner top edge
Line(6) = {6, 7}; // Right edge of top beam
Line(7) = {7, 8}; // GAMMA_D (Loading)
Line(8) = {8, 9}; 
Line(9) = {9, 1}; // LEFT_WALL (Clamped)

// Void boundaries (Full height right-side extension)
Line(10) = {2, 10}; // Void bottom right
Line(11) = {10, 11}; // VOID_RIGHT (Free edge)
Line(12) = {11, 7}; // Void top right

// -----------------------------------------------------------------------------
// Surfaces (Curve Loops)
// -----------------------------------------------------------------------------
// C-Shape Solid (Counter-Clockwise)
Curve Loop(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9}; 
Plane Surface(1) = {1}; 

// Void Domain (Counter-Clockwise, enveloping inner cavity and right edges of both beams)
Curve Loop(2) = {10, 11, 12, -6, -5, -4, -3, -2};
Plane Surface(2) = {2};

// -----------------------------------------------------------------------------
// Physical Groups Output
// -----------------------------------------------------------------------------
Coherence;

Physical Surface("SOLID", 100) = {1};
Physical Surface("VOID", 101) = {2};

Physical Curve("LEFT_WALL", 200) = {9};
Physical Curve("GAMMA_D", 201) = {7};
Physical Curve("VOID_RIGHT", 202) = {11};

// Attempt to merge into quad elements
Recombine Surface {1, 2};
