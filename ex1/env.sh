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

export MFEM_DIR="/home/th202608/codes/mfem/master/install/lib/cmake/mfem"
export MFEMMGIS_DIR="/home/th202608/codes/mfem-mgis/master/install/share/mfem-mgis/cmake"
export MFrontGenericInterface_DIR="/home/th202608/codes/mgis/master/install-python-3.7/share/mgis/cmake"

pathadd "/home/th202608/codes/tfel/master/install-python-3.7/bin"
pathadd ""
pathadd "/usr/bin"
export PATH

ldadd "/home/th202608/codes/tfel/master/install-python-3.7/lib"
ldadd "/home/th202608/codes/mgis/master/install-python-3.7/share/mgis/cmake/../../../lib"
export LD_LIBRARY_PATH

export CC=/usr/bin/cc
export CXX=/usr/bin/c++




