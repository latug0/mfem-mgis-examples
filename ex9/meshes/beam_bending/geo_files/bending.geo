// -----------------------------------------------------------------------------
// Parameters 
// -----------------------------------------------------------------------------
factor = 1.0;

// Poutre
length = 2.0 * factor;
width  = 0.1 * factor;
hpo    = 0.55 * factor; // Distance poutre - origine

// Premier corps (en haut a gauche)
dxB = 0.25 * length;
dyB = 6.0 * width;
hpb = 0.05 * length; // Distance poutre - premier corps

// Deuxieme corps (en bas a droite)
diam   = 0.3 * length;
hcb    = 2.0 * width;
dpb    = 1.0 * length; // Distance deuxieme corps - origine
radius = diam / 2.0;

// Mesh Density
If (!Exists(ref))
  ref = 1;
EndIf

d1 = (0.01 / ref) * length;
lc = d1;

// Logical Test
BOOL01 = (diam + dpb >= length);

// -----------------------------------------------------------------------------
// Points Construction
// -----------------------------------------------------------------------------

// Points principaux de la poutre
Point(1) = {0, hpo, 0, lc};             // POU00
Point(2) = {length, hpo, 0, lc};        // POU01
Point(3) = {length, hpo+width, 0, lc};  // POU02
Point(4) = {0, hpo+width, 0, lc};       // POU03
Point(6) = {dxB, hpo+width, 0, lc};     // POU05

If (BOOL01)
    Point(5) = {length-diam, hpo, 0, lc}; // POU04
Else
    Point(5) = {dpb, hpo, 0, lc};         // POU04
    Point(7) = {dpb+diam, hpo, 0, lc};    // POU06
EndIf

// Points pour le Corps 1
ybod = hpo + width + hpb + dyB;
pden = (dyB * dyB) - (dxB * dxB);
pcer = pden / (2.0 * dxB);

Point(201) = {0, ybod-dyB, 0, lc};      // PBOD101
Point(202) = {dxB, ybod, 0, lc};        // PBOD102
Point(203) = {0, ybod, 0, lc};          // PBOD103
Point(204) = {-pcer, ybod, 0, lc};      // PBOD104 (Center)

// Points pour le Corps 2
If (BOOL01)
    Point(101) = {length-diam, 0, 0, lc};      // PBOD201
    Point(102) = {length, 0, 0, lc};           // PBOD202
    Point(103) = {length, hcb, 0, lc};         // PBOD203
    Point(104) = {length-diam, hcb, 0, lc};    // PBOD204
    Point(105) = {length-radius, hcb, 0, lc};  // PBOD205 (Center)
Else
    Point(101) = {dpb, 0, 0, lc};              // PBOD201
    Point(102) = {dpb+diam, 0, 0, lc};         // PBOD202
    Point(103) = {dpb+diam, hcb, 0, lc};       // PBOD203
    Point(104) = {dpb, hcb, 0, lc};            // PBOD204
    Point(105) = {dpb+radius, hcb, 0, lc};     // PBOD205 (Center)
EndIf

// -----------------------------------------------------------------------------
// Curves & Arcs Construction
// -----------------------------------------------------------------------------

// Lignes de la Poutre
Line(1) = {1, 5}; // LPOUB
Line(3) = {2, 3}; // LPOUD
Line(4) = {3, 6}; // LPOUH
Line(5) = {6, 4}; // LPOUV
Line(6) = {4, 1}; // LPOUG

If (BOOL01)
    Line(2) = {5, 2}; // LPOUI
Else
    Line(20) = {5, 7}; // LPOUI
    Line(21) = {7, 2}; // LPOUK
EndIf

// Arcs pour les Corps
Circle(11) = {201, 204, 202}; // LBOD1D
Circle(10) = {103, 105, 104}; // LBOD2H (Traced right-to-left to ensure CCW upward bump)

// Boundary lines for Body 2 (exported but not strictly defining meshed domains)
Line(30) = {103, 102}; // LBOD2D
Line(31) = {102, 101}; // LBOD2G

// Lignes pour les Voids
Line(12) = {4, 201}; // LVOI1G (POU03 to PBOD101)
Line(13) = {202, 6}; // LVOI1D (PBOD102 to POU05)

If (BOOL01)
    Line(14) = {104, 5}; // LVOI2G (PBOD204 to POU04)
    Line(15) = {2, 103}; // LVOI2D (POU01 to PBOD203)
Else
    Line(14) = {104, 5}; // LVOI2G (PBOD204 to POU04)
    Line(15) = {7, 103}; // LVOI2D (POU06 to PBOD203)
EndIf

// -----------------------------------------------------------------------------
// Surfaces (Curve Loops must be continuous & correctly oriented)
// -----------------------------------------------------------------------------

// POUTRE
If (BOOL01)
    Curve Loop(1) = {1, 2, 3, 4, 5, 6}; 
Else
    Curve Loop(1) = {1, 20, 21, 3, 4, 5, 6};
EndIf
Plane Surface(1) = {1}; 

// VOID1
Curve Loop(2) = {12, 11, 13, 5};
Plane Surface(2) = {2};

// VOID2
If (BOOL01)
    Curve Loop(3) = {14, 2, 15, 10}; 
Else
    Curve Loop(3) = {14, 20, 15, 10};
EndIf
Plane Surface(3) = {3};

// -----------------------------------------------------------------------------
// Physical Groups Output (Equivalents of Cast3M domains & boundaries)
// -----------------------------------------------------------------------------

// CRITICAL: Mathematically fuse the domains
Coherence;

// Domains
Physical Surface("POUTRE", 100) = {1};
Physical Surface("VOID1", 101) = {2};
Physical Surface("VOID2", 102) = {3};

// Boundaries - POUTRE
Physical Curve("LPOUB", 201) = {1};
If (BOOL01)
    Physical Curve("LPOUI", 202) = {2};
Else
    Physical Curve("LPOUI", 202) = {20};
    Physical Curve("LPOUK", 203) = {21};
EndIf
Physical Curve("LPOUD", 204) = {3};
Physical Curve("LPOUH", 205) = {4};
Physical Curve("LPOUV", 206) = {5};
Physical Curve("LPOUG", 207) = {6};

// Boundaries - Voids
Physical Curve("LBOD1D", 208) = {11};
Physical Curve("LBOD2H", 209) = {10};
Physical Curve("LVOI1G", 212) = {12};
Physical Curve("LVOI1D", 213) = {13};
Physical Curve("LVOI2G", 214) = {14};
Physical Curve("LVOI2D", 215) = {15};

