# Installation Guide on Topaze/CCRT

This guide provides step-by-step instructions for setting up your environment on Topaze/CCRT and installing the necessary software. Follow these steps to get started.

## Create a new directory and useful paths

```
mkdir topaze-dir && cd topaze-dir
export MY_DIR=$PWD
export MY_LOG=YOURLOGIN
export MY_DEST=/ccc/scratch/cont002/den/YOURLOGIN/mini-test
```

## Download Spack

How to download Spack: 

```
cd $MY_DIR
git clone https://github.com/spack/spack.git
export SPACK_ROOT=$PWD/spack
```

Before proceeding, make sure to source Spack and clear your local ~/.spack repository (warning).

```
rm -r ~/.spack
source ${SPACK_ROOT}/share/spack/setup-env.sh
```

## Create a Spack Mirror on Your Machine (Local)

Firstly, you need to get the mmm spack repository.

```
git clone https://github.com/thelfer/mfem-mgis
spack repo add $PWD/mfem-mgis/spack_repo
```

Now, you will create a spack mirror and a boostrap directory.

```
spack bootstrap mirror --binary-packages my_bootstrap
spack mirror create -d mirror-mmm -D mmm+mpi+suite-sparse%gcc@11.1.0
```

It's possible that you will need some packages in your mirror, you can specify them with the following command:

```
spack mirror create -d mirror-mmm -D mmm+mpi+suite-sparse zlib ca-certificates-mozilla zlib-ng util-macros pkgconf findutils libpciaccess libedit libxcrypt bison libevent numactl
```

## Copy Data to Topaze

You'll need to copy the following files to Topaze:
- spack
  spack
- mfem-mgis
- mfem-mgis-example

Create an archive for these files:

```
cd $MY_DIR
tar cvf archive.tar.gz mfem-mgis/ mfem-mgis-examples/ mirror-mmm/ spack/ my_bootstrap/
scp archive.tar.gz $MY_LOG@topaze.ccc.cea.fr:$MY_DEST/
```

## Load Topaze modules
 
Log on `Topaze`: 

```
ssh -Y $MY_LOG@topaze.ccc.cea.fr
```

Load the required modules on Topaze:

```
module load gnu/11.1.0
module load mpi
```
## Install mfem-mgis on Topaze

Note that the installation is performed in your scratch directory, and files are automatically removed after 3 months.

### Setup spack

```
cd $MY_DEST
tar xvf archive.tar.gz
source $PWD/spack/share/spack/setup-env.sh
spack bootstrap reset -y
spack bootstrap add --scope=site --trust local-binaries $PWD/my_bootstrap/metadata/binaries/
spack bootstrap disable --scope=site github-actions-v0.5
spack bootstrap disable --scope=site github-actions-v0.4
spack bootstrap disable --scope=site spack-install
spack bootstrap root $PWD/spack/bootstrap
spack repo add mfem-mgis/spack_repo/
spack bootstrap now
spack bootstrap status
```

### Export SPACK Variables

To use MFront, you need to export some SPACK variables. Please execute the following commands:

```
export CC='gcc'
export CXX='g++'
export FC='mpifort'
export OMPI_CC='gcc'
export OMPI_CXX='g++'
export OMPI_FC='gfortran'
```

### Install MFEM-MGIS

```
spack repo add $PWD/mfem-mgis/spack_repo/
spack mirror add MMM $PWD/mirror-mmm/
```

Run installation:

```
module load gnu/11.1.0 mpi hwloc cmake
spack compiler find
spack external find hwloc
spack external find cmake
spack external find openssh
spack external find openmpi
spack install mmm+mpi+suite-sparse%gcc@11.1.0
```

### Install MFEM-MGIS-example on Topaze

Follow these steps to install mfem-mgis-example on Topaze:

```
cd mfem-mgis-example
mkdir build && cd build
spack load mmm
export MFEMMGIS_DIR=`spack location -i mmm`/share/mfem-mgis/cmake/
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
