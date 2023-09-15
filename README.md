# mfem-mgis-data

This repository provides user oriented material for the external `mfem-mgis` [repository](https://github.com/thelfer/mfem-mgis).
You will find here meshes and a collection of use cases for the `mfem-mgis` library.

## How to install

### Installing MFEM-MGIS with Spack

Firstly, get the mfem-mgis spack repository.

```
git clone https://github.com/thelfer/mfem-mgis
spack repo add $PWD/mfem-mgis/spack_repo
```

Secondly, install mfem-mgis

```
spack install mmm+mpi+suite-sparse
```

Thirdly, load mfem-mgis

```
spack load mfem-mgis
```

Finally, build your examples:

```
mkdir build
cd build
make -j 4
```

## Test case description

| Name | Description | Directory
|--|--|--|
| TensileTest | TODO | ex1 |
| Ssna303     | This tutorial deals with a 2D (plane strain) tensile test on a notched beam modeled by finite-strain plastic behavior. A tutorial describing this simulation is available at: https://thelfer.github.io/mfem-mgis/web/tutorial.html | ex2 |
| Inclusions  | This example models several inclusions in a periodic cube under imposed macroscopic strain. | ex3 |
| Ssna303_3d  | This tutorial deals with a 3D tensile test on a notched beam modeled by finite-strain plastic behavior. | ex4 |
| Satoh       | Modelling plate of length 1 in plane strain clamped on the left and right boundaries and submitted to a parabolic thermal gradient along the x-axis | ex5 |
| Rve-elastic | Simulation of a Representative Volume Element (RVE) with a non-linear elastic behavior law. A geometry mesh is provided : "inclusions_49.geo". The mesh can be generated using the following command: gmsh -3 `inclusions_49.geo`. By modifying the parameters within the `.geo` file, such as the number of spheres and the size of the element mesh, you can control and customize the simulation accordingly  | ex6 |
| Mox2        | Simulation of a Representative Volume Element (RVE) Mixed OxideFuels  with a viscoplastic behavior law. A mesh with one inclusion is provided : inclusion.msh. More information in ex7/README.md  | ex7 |

