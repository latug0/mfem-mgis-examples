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

export MFEM_DIR="/home/gl224549/Work/Dev/helix/mfem/install/lib/cmake/mfem"
export MFEMMGIS_DIR="/home/gl224549/Work/Dev/mfem-mgis/build_parallel/install/share/mfem-mgis/cmake"
export MFrontGenericInterface_DIR="/home/gl224549/spack/opt/spack/linux-debian9-skylake_avx512/gcc-9.2.0/mgis-master-dxbc5zylhv25ten4equum6z6kwt7sf7j/share/mgis/cmake"

pathadd "/home/gl224549/spack/opt/spack/linux-debian9-skylake_avx512/gcc-9.2.0/tfel-master-gudzydfnwltqnohtkihohehgsmrfe2jm/bin"
pathadd "/home/gl224549/spack/opt/spack/linux-debian9-skylake_avx512/gcc-9.2.0/openmpi-3.1.5-fa6b2jot7wirv5h44rom6lbxl2xue36c/bin"
pathadd "/home/gl224549/spack/opt/spack/linux-debian9-skylake_avx512/gcc-9.2.0/cmake-3.16.2-sda3osp22ts55fu4mqdv4og4s7m37t4y/bin"
export PATH

ldadd "/home/gl224549/spack/opt/spack/linux-debian9-skylake_avx512/gcc-9.2.0/tfel-master-gudzydfnwltqnohtkihohehgsmrfe2jm/lib"
ldadd "/home/gl224549/spack/opt/spack/linux-debian9-skylake_avx512/gcc-9.2.0/mgis-master-dxbc5zylhv25ten4equum6z6kwt7sf7j/share/mgis/cmake/../../../lib"
export LD_LIBRARY_PATH

export CC=/home/gl224549/spack/opt/spack/linux-debian9-skylake_avx512/gcc-9.2.0/gcc-9.2.0-3b2g3pktyfgjw7b4jmpeunbqoplys3ib/bin/gcc
export CXX=/home/gl224549/spack/opt/spack/linux-debian9-skylake_avx512/gcc-9.2.0/gcc-9.2.0-3b2g3pktyfgjw7b4jmpeunbqoplys3ib/bin/g++




