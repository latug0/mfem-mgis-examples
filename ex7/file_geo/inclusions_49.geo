// -----------------------------------------------------------------------------
//
//  Close to the Gmsh GEO tutorial 16
//
//  Constructive Solid Geometry, OpenCASCADE geometry kernel
//
// -----------------------------------------------------------------------------

// Please use GMSH version 4.7
SetFactory("OpenCASCADE");

nb = 50;
radius = 0.12;
default_edge_size = .015;
box_edge_size=1.;

// margin is taken a bit larger than radius to avoid the sphere beeing close to boundary
margin = .02;

//-------------------------------
Box(1) = {0,0,0, box_edge_size,box_edge_size,box_edge_size};
x = {};
y = {};
z = {};
r = {};
lc = {};
vol_sphere = 0.;
base = 2; // has to be strictly larger than 1
nbfails = 0;
t = 0;
For count In {0:1000*nb-1}
  If (t < nb)
    r[t]=radius/6+Rand(radius(1-1/6));
    x[t]=r[t]+Rand(1.-2*(r[t]+margin));
    y[t]=r[t]+Rand(1.-2*(r[t]+margin));
    z[t]=r[t]+Rand(1.-2*(r[t]+margin));
    intersect=0;
    If (x[t]-r[t] < margin || x[t]+r[t]+margin > box_edge_size || y[t]-r[t] < margin || y[t]+r[t]+margin > box_edge_size || z[t]-r[t] < margin || z[t]+r[t]+margin > box_edge_size)
      intersect=1;
    Else
      For i In {0:t-1}
        dist = Sqrt((x[t]-x[i])*(x[t]-x[i])+(y[t]-y[i])*(y[t]-y[i])+(z[t]-z[i])*(z[t]-z[i]));
        If (dist < r[t]+r[i])
          intersect=1;
        EndIf
      EndFor
    EndIf      
    If (intersect == 1)
      nbfails=nbfails+1;
    Else
      Sphere(base + t) = {x[t],y[t],z[t],r[t]};
      local_vol = 4.*Pi*r[t]*r[t]*r[t]/3.;
      vol_sphere = vol_sphere + local_vol;
      Printf("sphere t=%5g rad=%12g vol=%12g cumulated_vol=%g", t, r[t], local_vol,vol_sphere);
      // Here, we specifiy the edge size for this sphere
      lc[t] = r[t]/3.;
      t = t + 1;
    EndIf
  EndIf
EndFor
If (t != nb)
  Printf("Fail to create the precscribed number of spheres");
  Exit;
EndIf
// Delete all objects and create new ones by BooleanFragment operator
v() = BooleanFragments{ Volume{1}; Delete; }{ Volume{base : base+nb-1}; Delete; };

// Assign a mesh size to all the points of all the volumes:
MeshSize{ PointsOf{ Volume{:}; } } = default_edge_size;

For t In {0:nb-1}
  MeshSize{ PointsOf{ Volume{v(t)}; } } = lc[t];
EndFor
Printf("Theoretical volumic fraction %g",vol_sphere/(box_edge_size*box_edge_size*box_edge_size));
Physical Volume(1) = v(#v()-1);
Physical Volume(2) = {v(0):v(nb-1)};


eps = 1e-4;
// Select the corner point by searching for it geometrically:
sxmin() = Surface In BoundingBox{0.-eps, 0.-eps, 0.-eps,
                                 0+eps, 1.+eps, 1.+eps};
sxmax() = Surface In BoundingBox{1.-eps, 0.-eps, 0.-eps,
                                 1+eps, 1.+eps, 1.+eps};

Periodic Surface {sxmax()} = {sxmin()} Translate {1,0,0}; 								 

symin() = Surface In BoundingBox{0.-eps, 0.-eps, 0.-eps,
                                 1+eps, 0.+eps, 1.+eps};
symax() = Surface In BoundingBox{0.-eps, 1.-eps, 0.-eps,
                                 1+eps, 1.+eps, 1.+eps};
								 
Periodic Surface {symax()} = {symin()} Translate {0,1,0}; 		

szmin() = Surface In BoundingBox{0.-eps, 0.-eps, 0.-eps,
                                 1+eps, 1.+eps, 0.+eps};
szmax() = Surface In BoundingBox{0.-eps, 0.-eps, 1.-eps,
                                 1+eps, 1.+eps, 1.+eps};
								 
Periodic Surface {szmax()} = {szmin()} Translate {0,0,1}; 							 

Physical Surface (1) = sxmin();								 
Physical Surface (2) = sxmax();								 
Physical Surface (3) = symin();								 
Physical Surface (4) = symax();								 
Physical Surface (5) = szmin();								 
Physical Surface (6) = szmax();								 
								 

Mesh.MshFileVersion = 2.2;
