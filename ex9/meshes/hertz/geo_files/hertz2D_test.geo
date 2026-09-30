// -----------------------------------------------------------------------------
// Parameters & Geometry for 2D Hertz Contact (Cylinder on Flat)
// -----------------------------------------------------------------------------
R  = 1; // Indenter radius
W  = 2*R;  // Half-width of plane (using symmetry)
W = Max(W,R);
Hf = 0.9*R;  // Foundation/plane depth
void_gap = 0.05*R; // Gap between indenter and foundation
arc_ratio = 0.5; // Percentage of arc length covered by third Medium

lc = 0.05;  // Mesh characteristic length
n_contact = 40; // Number of nodes in contact area
N = n_contact - 1; // Number of elements in contact area

// Foundation
Point(1) = {0, -Hf, 0, lc};
Point(2) = {W, -Hf, 0, lc};
Point(3) = {W, 0, 0, lc};
Point(4) = {R*Sin(arc_ratio*Pi/2),0,0,lc}; 
Point(5) = {0, 0, 0, lc}; 

Line(1) = {1, 2}; // FOUND_BOT
Line(2) = {2, 3}; // FOUND_RIGHT
Line(3) = {3, 4}; // FOUND_TOP

// --- Align Line 4 nodes with Circle 6 to ensure vertical elements ---
// We generate exact node coordinates matching the projected arc and 
// chain them into a continuous curve loop list.
pt_list[0] = 4;
For i In {1 : N-1}
    // Angle uniformly decreases from max contact angle down to 0
    theta = (arc_ratio*Pi/2) * (1 - i/N);
    x = R*Sin(theta);
    Point(100+i) = {x, 0, 0, lc};
    pt_list[i] = 100+i;
EndFor
pt_list[N] = 5;

line4_list[] = {};
For i In {1 : N}
    Line(400+i) = {pt_list[i-1], pt_list[i]};
    line4_list[i-1] = 400+i;
    Transfinite Curve {400+i} = 2; // 2 nodes per segment (1 element)
EndFor
// --------------------------------------------------------------------

Line(5) = {5, 1}; // SYM_FOUND (x=0)

// Curve Loop 1 uses the array of segments instead of a single Line(4)
Curve Loop(1) = {1, 2, 3, line4_list[], 5};
Plane Surface(1) = {1}; 

// Void (Third Medium Gap)
Point(6) = {0,void_gap,0,lc};
Point(7) = {R*Sin(arc_ratio*Pi/2) , void_gap + R*(1-Cos(arc_ratio*Pi/2)), 0, lc};
Point(8) = {R, void_gap+R, 0, lc}; 
Point(9) = {0, R+void_gap, 0, lc}; 

Circle(6) = {6, 9, 7};    // IND_CONTACT (Arc) 
Circle(7) = {7, 9, 8};    // IND_BOT (Arc) 

Line(8) = {5,6}; // VOID_LEFT 
Line(9) = {7, 4}; // VOID_RIGHT

// Curve Loop 2 uses the array of segments to close the void
Curve Loop(2) = {8, 6, 9, line4_list[]}; 
Plane Surface(2) = {2};

// Indenter
Line(10) = {8, 9}; // IND_TOP
Line(11) = {9, 6}; // SYM_IND (x=0)

Curve Loop(3) = {6, 7, 10, 11}; 
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

// -----------------------------------------------------------------------------
// Quadrilateral Meshing Algorithms
// -----------------------------------------------------------------------------
Transfinite Curve {6} = n_contact; // Top contact arc
Transfinite Curve {8,9} = 2;       // 1 element width for void

// Because the bottom boundary is now composed of multiple segments, 
// the corners must be explicitly defined for the Transfinite Surface mapping to work.
Transfinite Surface {2} = {5, 6, 7, 4}; 
Recombine Surface {2};