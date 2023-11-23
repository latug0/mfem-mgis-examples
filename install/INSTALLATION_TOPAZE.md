# Installation Guide on Topaze/CCRT

This guide provides step-by-step instructions for setting up your environment on Topaze/CCRT and installing the necessary software. Follow these steps to get started.

## Create a Spack Mirror on Your Host Machine
Before proceeding, make sure to source Spack.


```
spack bootstrap mirror --binary-packages my_bootstrap
spack mirror create -d mirror-mmm -D mmm+mpi+suite-sparse
```

## Copy Data to Topaze

You'll need to copy the following files to Topaze:
- spack
- mfem-mgis
- mfem-mgis-example

Create an archive for these files:

```
tar cvf archive.tar.gz mfem-mgis/ mfem-mgis-examples/ mirror-mmm/ spack/ my_bootstrap/
scp archive.tar.gz yourlogin@topaze.ccc.cea.fr:.
```

## Load Topaze modules
 
Load the required modules on Topaze:

```
module load gnu/12.2.0
module load mpi/openmpi/4.0.5
```

## Install mfem-mgis on Topaze

Installation is performed in your scratch directory, and files are automatically removed after 3 months.

```
cd 
mv archive.tar.gz $CCCSCRATCHDIR/your-directory
cd $CCCSCRATCHDIR/your-directory
tar xvf archive.tar.gz
spack bootstrap add --trust local-sources my_bootstrap/metadata/sources/
spack bootstrap add --trust local-binaries my_bootstrap/metadata/binaries/
spack repo add mfem-mgis/spack_repo/
spack install mmm+mpi+suite-sparse
```

## Export SPACK Variables

To use MFront, you need to export some SPACK variables. Please execute the following commands:

```
export CC='gcc'
export CXX='g++'
export FC='mpifort'
export OMPI_CC='gcc'
export OMPI_CXX='g++'
export OMPI_FC='gfortran'
```

Don't worry; these variables are set up even though we won't use Fortran.

## Install mfem-mgis-example on Topaze

Follow these steps to install mfem-mgis-example on Topaze:

```
cd mfem-mgis-example
mkdir build && cd build
cmake ..
make -j 10
ctest
```

## How to run an example (ex8)

There are two ways to run an example, such as ex8:

### Using ccc_mprun

To run an example using ccc_mprun with 32 processes and 1 core per process, execute the following command:

```
ccc_mprun -n 32 -c 1 -pmilan ./uniaxial-elastic
```

### Using ccc_msub 

TODO

## Troubleshooting

### If you encounter Spack errors due to missing packages, consider the following possibilities:

Two possibilities:

1) Check if the package is already installed on Topaze by running:

```
spack external find your-package
```
If the package is found, you can use it directly.

2) If the package is not installed on Topaze, you can add its sources to your mirror directory. If you are using an SSHFS mount, you can complete your mirror by executing the following command on your host machine: 

```
spack mirror create -d your-mirror/ -D your-package
```

For more questions about spack, see the spack documentation.
