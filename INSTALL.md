---
title: Installation guide
author: Guillaume Latu, Thomas Helfer
date: 11/04/2021
lang: en-EN
link-citations: true
colorlinks: true
numbersections: true
toc: true
geometry:
  - margin=2cm
papersize: a4
figPrefixTemplate: "$$i$$"
tblPrefixTemplate: "$$i$$"
secPrefixTemplate: "$$i$$"
eqnPrefixTemplate: "($$i$$)"
---

This project uses [`cmake`](https://cmake.org/) as build system.

# Dependencies

List of prerequisites:

- [`MFEM`](https://mfem.org/)
- [`MGIS`](https://github.com/thelfer/MFrontGenericInterfaceSupport)
- [`MFEM-MGIS`](https://github.com/thelfer/mfem-mgis)

We refer to the `mfem-mgis` installation guide to install
prerequisites and `mfem-mgis` library.
In addition to that, you shall activate the MUMPS support or SUITESPARSE support in MFEM. Without this, you will not be able to activate all the examples.

Warning: to have access to all the examples, you should :

- go through the `mfem-mgis` install procedure and perform a `make install`.
- define an environment variable MFEMMGIS_DIR that contains the cmake configure file. For example with the command :
~~~~{.bash}
   export MFEMMGIS_DIR=${WORK_PLACE}/mfem-mgis/install/share/mfem-mgis/cmake
~~~~

# Check all exemples are correctly executed

~~~~{.bash}
$ mkdir build
$ cmake .. -DCMAKE_INSTALL_PREFIX=$PWD/install
$ make -j check && make install
~~~~

# Creating a simple example based on `mfem-mgis`

Through the `make install` command, some examples has been created in
your installation directory.

You can copy them elsewhere together with the `env.sh` file. The example
can be compiled either using the build systems `cmake` or`make`.

## Building the example using the `cmake` build-system

~~~~{.bash}
$ export INSTALLDIR=<your_mfemmgis_install_directory>
$ export TGDIR=<your_work_directory>
$ cd ${TGDIR}
$ cp -r ${INSTALLDIR}/share/mfem-mgis-examples/ex1 .
$ cp ${INSTALLDIR}/share/mfem-mgis-examples/env.sh ex1/
$ cd ex1
$ source env.sh
$ mkdir build
$ cd build
$ cmake ..
$ make
~~~~

The example may then be run as follows:

~~~~{.bash}
$ ./UniaxialTensileTest 
~~~~

You can then modify the source file and design your
own case of study.

## Building the example using the `make` build-system

~~~~{.bash}
$ export INSTALLDIR=<your_mfemmgis_install_directory>
$ export TGDIR=<your_work_directory>
$ cd ${TGDIR}
$ cp -r ${INSTALLDIR}/share/mfem-mgis-examples/ex1 .
$ cp ${INSTALLDIR}/share/mfem-mgis-examples/env.sh ex1/
$ cd ex1
$ source env.sh
$ make
~~~~

### Building in debug mode

The example and the `MFront` behaviour may be compiled in `debug` mode
by changing the call to make as follows:

~~~~{.bash}
$ make DEBUG=1
~~~~
