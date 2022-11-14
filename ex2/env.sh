#!/bin/bash
ldadd() {
    if [ -d "$1" ] && [[ ":$LD_LIBRARY_PATH:" != *":$1:"* ]]; then
        LD_LIBRARY_PATH="$1${LD_LIBRARY_PATH:+":$LD_LIBRARY_PATH"}"
    fi
}
pathadd() {
    if [ -d "$1" ] && [[ ":$PATH:" != *":$1:"* ]]; then
        PATH="$1${PATH:+":$PATH"}"
    fi
}

export MFEM_DIR="/ccc/cont002/home/den/pratraph/codes/mfem/install/lib/cmake/mfem"
export MFEMMGIS_DIR="/ccc/cont002/home/den/pratraph/codes/mfem-mgis/install/share/mfem-mgis/cmake"
export MFrontGenericInterface_DIR="/ccc/cont002/dsku/lautrec/home/user/den/pratraph/spack/opt/spack/linux-rhel8-zen2/gcc-11.1.0/mgis-master-taixbros6y6rkodhna7ns4kza2ductfn/share/mgis/cmake"

pathadd "/ccc/cont002/home/den/pratraph/spack/opt/spack/linux-rhel8-zen2/gcc-11.1.0/tfel-master-4e6ol7umoezxcvsprbdpuhyj5rna5y2f/bin"
pathadd "/ccc/products/openmpi-4.0.5/gcc--11.1.0/default/bin"
pathadd "/usr/bin"
export PATH

ldadd "/ccc/cont002/dsku/lautrec/home/user/den/pratraph/spack/opt/spack/linux-rhel8-zen2/gcc-11.1.0/tfel-master-4e6ol7umoezxcvsprbdpuhyj5rna5y2f/lib"
ldadd "/ccc/cont002/dsku/lautrec/home/user/den/pratraph/spack/opt/spack/linux-rhel8-zen2/gcc-11.1.0/mgis-master-taixbros6y6rkodhna7ns4kza2ductfn/share/mgis/cmake/../../../lib"
export LD_LIBRARY_PATH

export CC=/ccc/products/gcc-11.1.0/system/default/bin/gcc
export CXX=/ccc/products/gcc-11.1.0/system/default/bin/g++




