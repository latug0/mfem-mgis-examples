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
typeAlgo = sac_de_billes.AlgoRSA_3D() 
theSpheres = sac_de_billes.throwSpheres_3D( sac_de_billes.RSA, sac_de_billes.Tore, L, randomSeed, [[0.04,0.17]], [2], distMin)

#thePhases =  algo_3D.getPhases()
for sphere in theSpheres:
    sphere.phase = 2

sphInc = merope.SphereInclusions_3D()
sphInc.setLength(L)
sphInc.setSpheres(theSpheres)
#sphInc.setPhases(thePhases)

mi = merope.MultiInclusions_3D()
mi.setInclusions(sphInc)
mi.setMatrixPhase(1)


meshGenerator = merope.mesh.MeshGenerator()
meshGenerator.setMeshOrder(1)
meshGenerator.setMeshSize(0.01)
meshGenerator.setMultiInclusions(mi)
meshGenerator.write("634Spheres.geo")

