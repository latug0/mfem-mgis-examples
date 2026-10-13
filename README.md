# mfem-mgis-examples

This repository provides user oriented material for the external `mfem-mgis` [repository](https://github.com/thelfer/mfem-mgis).
You will find here meshes and a collection of use cases for the `mfem-mgis` library.

## How to install

### Installing MFEM-MGIS with Spack

```
git clone https://github.com/spack/spack.git
export SPACK_ROOT=$PWD/spack
source ${SPACK_ROOT}/share/spack/setup-env.sh
```

Firstly, install the development version of mfem-mgis, which the examples follow. It is
provided by the `develop` branch of the Spack packages repository, used by default by this
clone of Spack. With a release of Spack, first switch to this branch with
`spack repo update builtin --branch develop`.

```
spack install mfem-mgis@master
```

Secondly, load mfem-mgis

```
spack load mfem-mgis
```

Finally, build and run your examples:

```
cd mfem-mgis-examples
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 4
ctest --test-dir build
```

The tests compute the complete simulations, except in the Debug and Coverage builds, where they
only compute their beginning to be faster. The `MFEM_MGIS_EXAMPLES_TEST_MODE` option, set to `full`
or `restricted`, selects the test mode explicitly:

```
cmake -B build -DCMAKE_BUILD_TYPE=Release -DMFEM_MGIS_EXAMPLES_TEST_MODE=restricted
```

For an installation on a supercomputer without internet please follow the procedure described here for MFEM-MGIS installation: https://thelfer.github.io/mfem-mgis/installation_guide/installation_guide.html#installation-guide-on-topaze-ccrt-of-mfem-mgis-examples

## Test case description

| Name | Description | Directory
|--|--|--|
| TensileTest | Cyclic tension-compression of a unit cube made of a plastic material with linear isotropic hardening. | ex1 |
| Ssna303     | 2D (plane strain) tensile test on a notched beam with a finite-strain plastic behaviour, described in the [tutorial](https://thelfer.github.io/mfem-mgis/user_guide/tutorial.html) of mfem-mgis. | ex2 |
| TwoLayerCube | Periodic cube made of two elastic layers under an imposed macroscopic strain. | ex3 |
| Ssna303_3d  | 3D tensile test on a notched beam with a finite-strain plastic behaviour, solved with MUMPS, hypre or PETSc. | ex4 |
| Satoh       | Thermo-elastic plate in plane strain, clamped on its left and right boundaries and subjected to a parabolic temperature profile. | ex5 |
| Rve-elastic | Periodic Representative Volume Element (RVE) made of two materials with a Saint Venant-Kirchhoff hyperelastic behaviour under an imposed macroscopic deformation gradient: a two-layer cube, or spherical inclusions. | ex6 |
| Mox2        | Representative Volume Element (RVE) of a mixed oxide fuel: elastic inclusions in a viscoplastic matrix under an imposed macroscopic strain. More information in ex7/README.md | ex7 |
| RJH plate   | Thermomechanical simulation of one ring of a RJH fuel assembly. The fuel is U3Si2, the cladding and stiffeners are ALFENI. The heat transfer is non-linear with convective boundary conditions. The mechanics is in finite strain with thermal expansion, irradiation creep, swelling and plasticity. Both are strongly coupled at each time step. More information in ex8/README.md | ex8 |
