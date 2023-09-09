# -*- coding:utf8 -*-
#
# Author: M. Josien
# Date: 27/09/2021
#
# Copyright : see License.txt
#
# Build a .geo file to be meshed by gmsh
# The geometry is made of spheres

import sac_de_billes
import merope


L = [1, 1, 1]
distMin = 0.001
randomSeed = 0
#rad = 0.343653069
rad = 0.34365
density = 0.17
typeAlgo = sac_de_billes.AlgoRSA_3D() 
theSpheres = sac_de_billes.throwSpheres_3D( sac_de_billes.RSA, sac_de_billes.Tore, L, randomSeed, [[rad,density]], [2], distMin)
for sphere in theSpheres:
    sphere.phase = 2

sphInc = merope.SphereInclusions_3D()
sphInc.setLength(L)
sphInc.setSpheres(theSpheres)

mi = merope.MultiInclusions_3D()
mi.setInclusions(sphInc)
mi.setMatrixPhase(1)


meshGenerator = merope.mesh.MeshGenerator()
meshGenerator.setMeshOrder(2)
meshGenerator.setMeshSize(0.5)
meshGenerator.setMultiInclusions(mi)
meshGenerator.write("OneSphere.geo")
